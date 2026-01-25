/*
 * Metal GPU Compute Shader for 8n+3 Prime Search
 *
 * Implements the FJ64_262K primality test on GPU using software-emulated
 * 128-bit arithmetic (Metal lacks native 128-bit integers).
 *
 * Key features:
 * - One threadgroup per n value, 256 threads testing different a values
 * - Largest-a-first strategy (smallest prime candidates first)
 * - Atomic early termination when solution found
 * - Trial division filters ~80% of composites before Miller-Rabin
 * - FJ64 hash table for deterministic 2-witness primality test
 */

#include <metal_stdlib>
using namespace metal;

/* ========================================================================== */
/* Software 128-bit Arithmetic                                                 */
/* ========================================================================== */

struct uint128_t {
    uint64_t lo;
    uint64_t hi;
};

/* 64x64 -> 128-bit multiplication using 4 32x32 multiplies */
inline uint128_t mul64x64(uint64_t a, uint64_t b) {
    uint64_t a_lo = a & 0xFFFFFFFF;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = b & 0xFFFFFFFF;
    uint64_t b_hi = b >> 32;

    uint64_t p0 = a_lo * b_lo;
    uint64_t p1 = a_lo * b_hi;
    uint64_t p2 = a_hi * b_lo;
    uint64_t p3 = a_hi * b_hi;

    uint64_t mid = (p0 >> 32) + (p1 & 0xFFFFFFFF) + (p2 & 0xFFFFFFFF);
    uint64_t carry1 = mid >> 32;

    return {(p0 & 0xFFFFFFFF) | (mid << 32), p3 + (p1 >> 32) + (p2 >> 32) + carry1};
}

/* 128-bit addition with carry */
inline uint128_t add128(uint128_t a, uint128_t b) {
    uint64_t lo = a.lo + b.lo;
    uint64_t carry = (lo < a.lo) ? 1ULL : 0ULL;
    return {lo, a.hi + b.hi + carry};
}

/* ========================================================================== */
/* Montgomery Arithmetic                                                       */
/* ========================================================================== */

/* Compute n' = -n^(-1) mod 2^64 using Newton's method */
inline uint64_t montgomery_inverse(uint64_t n) {
    uint64_t x = n;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    return (uint64_t)(-(int64_t)x);
}

/* Montgomery reduction: compute t * r^(-1) mod n where r = 2^64 */
inline uint64_t montgomery_reduce(uint128_t t, uint64_t n, uint64_t n_inv) {
    uint64_t m = t.lo * n_inv;
    uint128_t mn = mul64x64(m, n);
    uint128_t sum = add128(t, mn);
    uint64_t u = sum.hi;
    return (u >= n) ? (u - n) : u;
}

/* Montgomery multiplication: a * b * r^(-1) mod n */
inline uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t n, uint64_t n_inv) {
    return montgomery_reduce(mul64x64(a, b), n, n_inv);
}

/* Compute r^2 mod n where r = 2^64 */
inline uint64_t montgomery_r_squared(uint64_t n) {
    uint64_t r = (1ULL << 63) % n;
    r = (r << 1);
    if (r >= n) r -= n;

    uint128_t r_sq = mul64x64(r, r);

    if (r_sq.hi == 0) {
        return r_sq.lo % n;
    }

    uint64_t hi_mod = r_sq.hi % n;
    uint128_t term = mul64x64(hi_mod, r);

    uint64_t lo_mod = r_sq.lo % n;
    uint128_t sum = {lo_mod, 0};
    sum = add128(sum, term);

    if (sum.hi > 0) {
        uint64_t sum_hi_mod = sum.hi % n;
        uint128_t term2 = mul64x64(sum_hi_mod, r);
        uint64_t sum_lo_mod = sum.lo % n;
        return (sum_lo_mod + (term2.lo % n)) % n;
    }

    return sum.lo % n;
}

/* Count trailing zeros in a 64-bit integer */
inline int ctz64(uint64_t x) {
    if (x == 0) return 64;
    int n = 0;
    if ((x & 0xFFFFFFFF) == 0) { n += 32; x >>= 32; }
    if ((x & 0x0000FFFF) == 0) { n += 16; x >>= 16; }
    if ((x & 0x000000FF) == 0) { n += 8;  x >>= 8; }
    if ((x & 0x0000000F) == 0) { n += 4;  x >>= 4; }
    if ((x & 0x00000003) == 0) { n += 2;  x >>= 2; }
    if ((x & 0x00000001) == 0) { n += 1; }
    return n;
}

/* ========================================================================== */
/* Miller-Rabin Witness Test                                                   */
/* ========================================================================== */

/* Miller-Rabin witness test using Montgomery multiplication */
inline bool mr_witness_montgomery(uint64_t n, uint64_t a, uint64_t n_inv, uint64_t r_sq) {
    if (a >= n) a = a % n;
    if (a == 0) return true;

    uint64_t d = n - 1;
    int r = ctz64(d);
    d >>= r;

    uint64_t one_m = montgomery_reduce({r_sq, 0}, n, n_inv);
    uint64_t neg_one_m = n - one_m;

    uint128_t a_rsq = mul64x64(a, r_sq);
    uint64_t a_m = montgomery_reduce(a_rsq, n, n_inv);

    uint64_t x_m = one_m;
    uint64_t base_m = a_m;
    uint64_t exp = d;

    while (exp > 0) {
        if (exp & 1) {
            x_m = montgomery_mul(x_m, base_m, n, n_inv);
        }
        base_m = montgomery_mul(base_m, base_m, n, n_inv);
        exp >>= 1;
    }

    if (x_m == one_m || x_m == neg_one_m)
        return true;

    for (int i = 1; i < r; i++) {
        x_m = montgomery_mul(x_m, x_m, n, n_inv);
        if (x_m == neg_one_m)
            return true;
        if (x_m == one_m)
            return false;
    }
    return false;
}

/* ========================================================================== */
/* FJ64 Hash and Primality Test                                                */
/* ========================================================================== */

/* FJ64 hash function - maps n to a bucket in [0, 262143] */
inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return (uint32_t)(x & 262143);
}

/* FJ64_262K primality test: exactly 2 Miller-Rabin tests */
inline bool is_prime_fj64(uint64_t n, device const uint16_t* fj64_bases) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);

    if (!mr_witness_montgomery(n, 2, n_inv, r_sq))
        return false;

    uint32_t hash_idx = fj64_hash(n);
    uint64_t witness2 = fj64_bases[hash_idx];

    return mr_witness_montgomery(n, witness2, n_inv, r_sq);
}

/* ========================================================================== */
/* Trial Division                                                              */
/* ========================================================================== */

constant uint32_t TRIAL_PRIMES[30] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};

/* Trial division: 0 = composite, 1 = small prime, 2 = needs MR */
inline int trial_division_check(uint64_t candidate) {
    /* Inline check first 7 primes */
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    if (candidate % 11 == 0) return (candidate == 11) ? 1 : 0;
    if (candidate % 13 == 0) return (candidate == 13) ? 1 : 0;
    if (candidate % 17 == 0) return (candidate == 17) ? 1 : 0;
    if (candidate % 19 == 0) return (candidate == 19) ? 1 : 0;

    /* Check remaining primes */
    for (int i = 7; i < 30; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* Full primality test with trial division */
inline bool is_candidate_prime(uint64_t candidate, device const uint16_t* fj64_bases) {
    int td = trial_division_check(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64(candidate, fj64_bases);
}

/* ========================================================================== */
/* Integer Square Root                                                         */
/* ========================================================================== */

inline uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    if (n == 1) return 1;

    int shift = (64 - clz(n)) >> 1;
    uint64_t x = 1ULL << (shift + 1);

    for (int i = 0; i < 8; i++) {
        uint64_t x_new = (x + n / x) >> 1;
        if (x_new >= x) break;
        x = x_new;
    }

    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;

    return x;
}

/* ========================================================================== */
/* Result Structure                                                            */
/* ========================================================================== */

struct SearchResult {
    uint64_t n;
    uint64_t a;
    uint64_t p;
    uint32_t found;
    uint32_t _pad;
};

/* ========================================================================== */
/* Search Kernel                                                               */
/* ========================================================================== */

kernel void search_kernel(
    device const uint64_t* n_values [[buffer(0)]],
    device const uint16_t* fj64_bases [[buffer(1)]],
    device SearchResult* results [[buffer(2)]],
    constant uint32_t& batch_size [[buffer(3)]],
    uint tg_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_index_in_threadgroup]],
    uint tg_size [[threads_per_threadgroup]]
) {
    if (tg_idx >= batch_size) return;

    uint64_t n = n_values[tg_idx];
    uint64_t N = 8 * n + 3;

    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    while (a_max >= 1 && a_max * a_max > N - 4) {
        a_max -= 2;
    }

    threadgroup atomic_uint found_flag;
    threadgroup uint64_t solution_a;
    threadgroup uint64_t solution_p;

    if (tid == 0) {
        atomic_store_explicit(&found_flag, 0, memory_order_relaxed);
        solution_a = 0;
        solution_p = 0;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint64_t num_a = (a_max + 1) / 2;

    for (uint64_t i = tid; i < num_a; i += tg_size) {
        if (atomic_load_explicit(&found_flag, memory_order_relaxed)) {
            break;
        }

        uint64_t a = a_max - 2 * i;
        if (a < 1) break;

        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;

        uint64_t candidate = (N - a_sq) >> 1;

        if (candidate >= 2 && is_candidate_prime(candidate, fj64_bases)) {
            uint expected = 0;
            if (atomic_compare_exchange_weak_explicit(&found_flag, &expected, 1,
                                                       memory_order_relaxed,
                                                       memory_order_relaxed)) {
                solution_a = a;
                solution_p = candidate;
            }
            break;
        }
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (tid == 0) {
        results[tg_idx].n = n;
        if (atomic_load_explicit(&found_flag, memory_order_relaxed)) {
            results[tg_idx].a = solution_a;
            results[tg_idx].p = solution_p;
            results[tg_idx].found = 1;
        } else {
            results[tg_idx].a = 0;
            results[tg_idx].p = 0;
            results[tg_idx].found = 0;
        }
    }
}

/* Alias for compatibility */
kernel void search_kernel_optimized(
    device const uint64_t* n_values [[buffer(0)]],
    device const uint16_t* fj64_bases [[buffer(1)]],
    device SearchResult* results [[buffer(2)]],
    constant uint32_t& batch_size [[buffer(3)]],
    uint tg_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_index_in_threadgroup]],
    uint tg_size [[threads_per_threadgroup]]
) {
    if (tg_idx >= batch_size) return;

    uint64_t n = n_values[tg_idx];
    uint64_t N = 8 * n + 3;

    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    while (a_max >= 1 && a_max * a_max > N - 4) {
        a_max -= 2;
    }

    threadgroup atomic_uint found_flag;
    threadgroup uint64_t solution_a;
    threadgroup uint64_t solution_p;

    if (tid == 0) {
        atomic_store_explicit(&found_flag, 0, memory_order_relaxed);
        solution_a = 0;
        solution_p = 0;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    uint64_t num_a = (a_max + 1) / 2;

    for (uint64_t i = tid; i < num_a; i += tg_size) {
        if (atomic_load_explicit(&found_flag, memory_order_relaxed)) {
            break;
        }

        uint64_t a = a_max - 2 * i;
        if (a < 1) break;

        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;

        uint64_t candidate = (N - a_sq) >> 1;

        if (candidate >= 2 && is_candidate_prime(candidate, fj64_bases)) {
            uint expected = 0;
            if (atomic_compare_exchange_weak_explicit(&found_flag, &expected, 1,
                                                       memory_order_relaxed,
                                                       memory_order_relaxed)) {
                solution_a = a;
                solution_p = candidate;
            }
            break;
        }
    }

    threadgroup_barrier(mem_flags::mem_threadgroup);

    if (tid == 0) {
        results[tg_idx].n = n;
        if (atomic_load_explicit(&found_flag, memory_order_relaxed)) {
            results[tg_idx].a = solution_a;
            results[tg_idx].p = solution_p;
            results[tg_idx].found = 1;
        } else {
            results[tg_idx].a = 0;
            results[tg_idx].p = 0;
            results[tg_idx].found = 0;
        }
    }
}
