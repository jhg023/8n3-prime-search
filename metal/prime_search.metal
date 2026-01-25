/*
 * Metal GPU Compute Shader for 8n+3 Prime Search
 *
 * Implements the FJ64_262K primality test on GPU using software-emulated
 * 128-bit arithmetic (Metal lacks native 128-bit integers).
 *
 * Key optimizations:
 * - One threadgroup per n value, 256 threads testing different a values
 * - Largest-a-first strategy (smallest prime candidates first)
 * - Atomic early termination when solution found
 * - Trial division filters ~80% of composites before Miller-Rabin
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

    /* Four 32x32 partial products */
    uint64_t p0 = a_lo * b_lo;  /* bits 0-63 */
    uint64_t p1 = a_lo * b_hi;  /* bits 32-95 */
    uint64_t p2 = a_hi * b_lo;  /* bits 32-95 */
    uint64_t p3 = a_hi * b_hi;  /* bits 64-127 */

    /* Combine middle terms with carry propagation */
    uint64_t mid = (p0 >> 32) + (p1 & 0xFFFFFFFF) + (p2 & 0xFFFFFFFF);
    uint64_t carry1 = mid >> 32;

    uint64_t hi = p3 + (p1 >> 32) + (p2 >> 32) + carry1;
    uint64_t lo = (p0 & 0xFFFFFFFF) | (mid << 32);

    return {lo, hi};
}

/* 128-bit addition with carry */
inline uint128_t add128(uint128_t a, uint128_t b) {
    uint64_t lo = a.lo + b.lo;
    uint64_t carry = (lo < a.lo) ? 1ULL : 0ULL;
    uint64_t hi = a.hi + b.hi + carry;
    return {lo, hi};
}

/* 128-bit subtraction (assumes a >= b) */
inline uint128_t sub128(uint128_t a, uint128_t b) {
    uint64_t borrow = (a.lo < b.lo) ? 1ULL : 0ULL;
    uint64_t lo = a.lo - b.lo;
    uint64_t hi = a.hi - b.hi - borrow;
    return {lo, hi};
}

/* Compare 128-bit values: returns -1 if a < b, 0 if a == b, 1 if a > b */
inline int cmp128(uint128_t a, uint128_t b) {
    if (a.hi < b.hi) return -1;
    if (a.hi > b.hi) return 1;
    if (a.lo < b.lo) return -1;
    if (a.lo > b.lo) return 1;
    return 0;
}

/* Check if 128-bit value is >= 64-bit value */
inline bool ge128_64(uint128_t a, uint64_t b) {
    return (a.hi > 0) || (a.lo >= b);
}

/* Right shift 128-bit value by 64 bits (returns high 64 bits) */
inline uint64_t shr128_64(uint128_t a) {
    return a.hi;
}

/* ========================================================================== */
/* Montgomery Arithmetic (Software 128-bit version)                            */
/* ========================================================================== */

/*
 * Compute n' = -n^(-1) mod 2^64 using Newton's method
 * Only valid for odd n (which is true for all prime candidates)
 */
inline uint64_t montgomery_inverse(uint64_t n) {
    uint64_t x = n;
    x *= 2 - n * x;  /* 4 bits */
    x *= 2 - n * x;  /* 8 bits */
    x *= 2 - n * x;  /* 16 bits */
    x *= 2 - n * x;  /* 32 bits */
    x *= 2 - n * x;  /* 64 bits */
    return (uint64_t)(-(int64_t)x);
}

/*
 * Montgomery reduction: compute t * r^(-1) mod n where r = 2^64
 * Input: t (128-bit), n (64-bit), n_inv = -n^(-1) mod 2^64
 * Output: t * 2^(-64) mod n
 */
inline uint64_t montgomery_reduce(uint128_t t, uint64_t n, uint64_t n_inv) {
    /* m = (t mod r) * n' mod r = t.lo * n_inv (low 64 bits only) */
    uint64_t m = t.lo * n_inv;

    /* mn = m * n (128-bit product) */
    uint128_t mn = mul64x64(m, n);

    /* u = (t + m*n) / r = (t + mn) >> 64 */
    uint128_t sum = add128(t, mn);
    uint64_t u = sum.hi;

    /* Final reduction: if u >= n, return u - n, else return u */
    return (u >= n) ? (u - n) : u;
}

/* Montgomery multiplication: a * b * r^(-1) mod n */
inline uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t n, uint64_t n_inv) {
    uint128_t product = mul64x64(a, b);
    return montgomery_reduce(product, n, n_inv);
}

/*
 * Compute r^2 mod n where r = 2^64
 * Uses the identity: r^2 mod n = (2^64 mod n)^2 mod n
 *
 * Since we can't compute 2^64 mod n directly (it overflows), we use:
 * 2^64 mod n = (2^63 mod n) * 2 mod n, adjusting if >= n
 * Then square the result and reduce.
 */
inline uint64_t montgomery_r_squared(uint64_t n) {
    /* Compute 2^64 mod n */
    /* 2^64 = (2^63) * 2, and 2^63 fits in uint64_t */
    uint64_t r = (1ULL << 63) % n;  /* 2^63 mod n */
    r = (r << 1);                    /* 2^64 mod n (may overflow if n is small) */
    if (r >= n) r -= n;

    /* Now compute r^2 mod n using 128-bit arithmetic */
    uint128_t r_sq = mul64x64(r, r);

    /* Reduce r_sq mod n using repeated subtraction (r_sq < n^2) */
    /* For n < 2^63, r_sq < 2^126, which is at most n * 2^63 */
    /* We need proper modular reduction here */

    /* Division by repeated subtraction is too slow, use a different approach */
    /* Compute (r mod n)^2 using schoolbook division */

    /* Since r < n (after reduction above), r^2 < n^2 */
    /* We need to compute r_sq mod n where r_sq is 128-bit and n is 64-bit */

    /* Barrett-style reduction: approximate division */
    /* For simplicity, use the fact that r_sq.hi < n (since r < n implies r^2 < n^2 < 2^64 * n for n > 2^32) */

    /* Actually, for n > 2^32, r < n implies r^2 can be up to ~2^126 */
    /* We need iterative reduction */

    /* Simplified approach: reduce high bits first */
    if (r_sq.hi == 0) {
        return r_sq.lo % n;
    }

    /* For r_sq.hi > 0, we have r_sq = r_sq.hi * 2^64 + r_sq.lo */
    /* r_sq mod n = ((r_sq.hi mod n) * (2^64 mod n) + r_sq.lo) mod n */
    /* where 2^64 mod n = r (computed above) */

    uint64_t hi_mod = r_sq.hi % n;
    uint128_t term = mul64x64(hi_mod, r);  /* hi_mod * (2^64 mod n) */

    /* Add lo part */
    uint64_t lo_mod = r_sq.lo % n;
    uint128_t sum = {lo_mod, 0};
    sum = add128(sum, term);

    /* Reduce again if needed */
    if (sum.hi > 0) {
        uint64_t sum_hi_mod = sum.hi % n;
        uint128_t term2 = mul64x64(sum_hi_mod, r);
        uint64_t sum_lo_mod = sum.lo % n;
        uint64_t result = (sum_lo_mod + (term2.lo % n)) % n;
        return result;
    }

    return sum.lo % n;
}

/* ========================================================================== */
/* Software CTZ (Count Trailing Zeros)                                         */
/* ========================================================================== */

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
/* Integer Square Root                                                         */
/* ========================================================================== */

/* Integer square root using Newton's method */
inline uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    if (n == 1) return 1;

    /* Initial estimate using leading zeros */
    int shift = (64 - clz(n)) >> 1;
    uint64_t x = 1ULL << (shift + 1);

    /* Newton's method iterations */
    for (int i = 0; i < 10; i++) {
        uint64_t x_new = (x + n / x) >> 1;
        if (x_new >= x) break;
        x = x_new;
    }

    /* Correct for any remaining error */
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;

    return x;
}

/* ========================================================================== */
/* Miller-Rabin Witness Test                                                   */
/* ========================================================================== */

/*
 * Miller-Rabin witness test using Montgomery multiplication
 * Returns true if n is a strong probable prime to base a
 */
inline bool mr_witness_montgomery(uint64_t n, uint64_t a, uint64_t n_inv, uint64_t r_sq) {
    if (a >= n) a = a % n;
    if (a == 0) return true;

    /* Compute d and r where n - 1 = d * 2^r */
    uint64_t d = n - 1;
    int r = ctz64(d);
    d >>= r;

    /* Montgomery form constants */
    /* one_m = 1 in Montgomery form = r mod n = reduce(r_sq) */
    uint64_t one_m = montgomery_reduce({r_sq, 0}, n, n_inv);
    uint64_t neg_one_m = n - one_m;  /* -1 in Montgomery form */

    /* Convert a to Montgomery form: a_m = a * r mod n = reduce(a * r_sq) */
    uint128_t a_rsq = mul64x64(a, r_sq);
    uint64_t a_m = montgomery_reduce(a_rsq, n, n_inv);

    /* Compute a^d mod n using Montgomery exponentiation */
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

    /* Check if x == 1 or x == n-1 (in Montgomery form) */
    if (x_m == one_m || x_m == neg_one_m) {
        return true;
    }

    /* Square r-1 times and check for n-1 */
    for (int i = 1; i < r; i++) {
        x_m = montgomery_mul(x_m, x_m, n, n_inv);
        if (x_m == neg_one_m) {
            return true;
        }
        if (x_m == one_m) {
            return false;  /* Composite */
        }
    }

    return false;  /* Composite */
}

/* ========================================================================== */
/* FJ64 Hash and Primality Test                                                */
/* ========================================================================== */

/* FJ64 hash function - maps n to a bucket in [0, 262143] */
inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return (uint32_t)(x & 262143);  /* & (2^18 - 1) */
}

/*
 * FJ64_262K primality test using Montgomery multiplication
 * Exactly 2 Miller-Rabin tests with optimized modular arithmetic
 * Assumes: n > 127, n is odd, n passed trial division
 */
inline bool is_prime_fj64(uint64_t n, device const uint16_t* fj64_bases) {
    /* Pre-compute Montgomery constants */
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);

    /* First witness: always 2 */
    if (!mr_witness_montgomery(n, 2, n_inv, r_sq)) {
        return false;
    }

    /* Second witness: from hash table */
    uint32_t hash_idx = fj64_hash(n);
    uint64_t witness2 = fj64_bases[hash_idx];

    return mr_witness_montgomery(n, witness2, n_inv, r_sq);
}

/* ========================================================================== */
/* Trial Division                                                              */
/* ========================================================================== */

/* 30 trial division primes (same as CPU version) */
constant uint32_t TRIAL_PRIMES[30] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};

/*
 * Trial division check
 * Returns: 0 = composite, 1 = is small prime, 2 = needs Miller-Rabin
 */
inline int trial_division_check(uint64_t candidate) {
    /* Check first 7 primes inline */
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    if (candidate % 11 == 0) return (candidate == 11) ? 1 : 0;
    if (candidate % 13 == 0) return (candidate == 13) ? 1 : 0;
    if (candidate % 17 == 0) return (candidate == 17) ? 1 : 0;
    if (candidate % 19 == 0) return (candidate == 19) ? 1 : 0;

    /* Check remaining primes (23-127) */
    for (int i = 7; i < 30; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }

    return 2;  /* Needs Miller-Rabin */
}

/* Full primality test */
inline bool is_candidate_prime(uint64_t candidate, device const uint16_t* fj64_bases) {
    int td = trial_division_check(candidate);
    if (td == 0) return false;  /* Composite */
    if (td == 1) return true;   /* Small prime */
    if (candidate <= 127) return true;  /* Small primes not caught by trial div */
    return is_prime_fj64(candidate, fj64_bases);
}

/* ========================================================================== */
/* Result Structure                                                            */
/* ========================================================================== */

/* Result for a single n value */
struct SearchResult {
    uint64_t n;          /* The n value */
    uint64_t a;          /* Solution a (0 if counterexample) */
    uint64_t p;          /* Solution p (0 if counterexample) */
    uint32_t found;      /* 1 if solution found, 0 if counterexample */
    uint32_t _pad;       /* Padding for alignment */
};

/* ========================================================================== */
/* Main Search Kernel                                                          */
/* ========================================================================== */

/*
 * Search kernel: one threadgroup per n value
 *
 * Each thread tests a different 'a' value (largest-first strategy).
 * First thread to find a valid (a, p) pair claims the result via atomic flag.
 *
 * Parameters:
 *   n_values: Array of n values to search
 *   fj64_bases: The 512KB hash table for FJ64 witness selection
 *   results: Output array of SearchResult structs
 *   batch_size: Number of n values in this batch
 */
kernel void search_kernel(
    device const uint64_t* n_values [[buffer(0)]],
    device const uint16_t* fj64_bases [[buffer(1)]],
    device SearchResult* results [[buffer(2)]],
    constant uint32_t& batch_size [[buffer(3)]],
    uint tg_idx [[threadgroup_position_in_grid]],
    uint tid [[thread_index_in_threadgroup]],
    uint tg_size [[threads_per_threadgroup]]
) {
    /* Check bounds */
    if (tg_idx >= batch_size) return;

    /* Get our n value */
    uint64_t n = n_values[tg_idx];
    uint64_t N = 8 * n + 3;

    /* Compute a_max (largest odd a such that a^2 < N) */
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;  /* Make odd */

    /* Ensure a_max^2 + 4 <= N (so p >= 2) */
    while (a_max >= 1 && a_max * a_max > N - 4) {
        a_max -= 2;
    }

    /* Shared memory for coordination */
    threadgroup atomic_uint found_flag;
    threadgroup uint64_t solution_a;
    threadgroup uint64_t solution_p;

    /* Thread 0 initializes shared state */
    if (tid == 0) {
        atomic_store_explicit(&found_flag, 0, memory_order_relaxed);
        solution_a = 0;
        solution_p = 0;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    /* Number of a values to test (a_max, a_max-2, ..., 1) */
    uint64_t num_a = (a_max + 1) / 2;

    /* Each thread tests a subset of a values */
    /* Thread tid tests a values: a_max - 2*tid, a_max - 2*(tid + tg_size), ... */
    for (uint64_t i = tid; i < num_a; i += tg_size) {
        /* Check if another thread already found a solution */
        if (atomic_load_explicit(&found_flag, memory_order_relaxed)) {
            break;
        }

        /* Compute a value (largest first) */
        uint64_t a = a_max - 2 * i;
        if (a < 1) break;

        /* Compute candidate p = (N - a^2) / 2 */
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;  /* p would be < 2 */

        uint64_t candidate = (N - a_sq) >> 1;

        /* Test primality */
        if (candidate >= 2 && is_candidate_prime(candidate, fj64_bases)) {
            /* Try to claim the solution */
            uint expected = 0;
            if (atomic_compare_exchange_weak_explicit(&found_flag, &expected, 1,
                                                       memory_order_relaxed,
                                                       memory_order_relaxed)) {
                /* We claimed it! Store the solution */
                solution_a = a;
                solution_p = candidate;
            }
            break;
        }
    }

    /* Wait for all threads */
    threadgroup_barrier(mem_flags::mem_threadgroup);

    /* Thread 0 writes the result */
    if (tid == 0) {
        results[tg_idx].n = n;
        if (atomic_load_explicit(&found_flag, memory_order_relaxed)) {
            results[tg_idx].a = solution_a;
            results[tg_idx].p = solution_p;
            results[tg_idx].found = 1;
        } else {
            results[tg_idx].a = 0;
            results[tg_idx].p = 0;
            results[tg_idx].found = 0;  /* Potential counterexample! */
        }
    }
}

/* ========================================================================== */
/* Verification Kernel (for testing 128-bit math)                              */
/* ========================================================================== */

/*
 * Test kernel for verifying 128-bit arithmetic correctness
 * Used during development/testing only
 */
kernel void test_mul64x64(
    device const uint64_t* a_values [[buffer(0)]],
    device const uint64_t* b_values [[buffer(1)]],
    device uint64_t* lo_results [[buffer(2)]],
    device uint64_t* hi_results [[buffer(3)]],
    uint idx [[thread_position_in_grid]]
) {
    uint128_t result = mul64x64(a_values[idx], b_values[idx]);
    lo_results[idx] = result.lo;
    hi_results[idx] = result.hi;
}
