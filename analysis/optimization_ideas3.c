/*
 * Optimization Ideas Benchmark - Part 3
 *
 * Test MR-related optimizations:
 * 1. Skip base-2 test when hash witness == 2
 * 2. Pre-compute Montgomery constants for small candidates
 * 3. Reduced r squarings in MR loop
 */

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/fj64_table.h"

/* Trial primes */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_TRIAL_PRIMES 30

#define ITERATIONS 10000000ULL

/* ========================================================================== */
/* Current Implementation (baseline)                                          */
/* ========================================================================== */

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline int trial_division_current(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_prime_fj64_current(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

static inline bool is_candidate_prime_current(uint64_t candidate) {
    int td = trial_division_current(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_current(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_current(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Optimization: Skip base-2 if hash witness == 2                             */
/* ========================================================================== */

static inline bool is_prime_fj64_skip2(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);

    uint32_t hash_witness = fj64_bases[fj64_hash(n)];

    /* If hash witness is 2, only do one test */
    if (hash_witness == 2) {
        return mr_witness_montgomery_cached(n, 2, n_inv, r_sq);
    }

    /* Otherwise, do both tests */
    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, hash_witness, n_inv, r_sq);
}

static inline bool is_candidate_prime_skip2(uint64_t candidate) {
    int td = trial_division_current(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_skip2(candidate);
}

static inline uint64_t find_solution_skip2(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_skip2(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Optimization: Check first 3 primes THEN compute hash                       */
/* ========================================================================== */

/* Idea: 3,5,7 catch ~50% of composites. Defer hash computation until needed */
static inline bool is_candidate_prime_defer_hash(uint64_t candidate) {
    /* Fast filter: 3, 5, 7 */
    if (candidate % 3 == 0) return candidate == 3;
    if (candidate % 5 == 0) return candidate == 5;
    if (candidate % 7 == 0) return candidate == 7;

    /* Check rest of trial primes */
    for (int i = 3; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return candidate == TRIAL_PRIMES[i];
        }
    }

    if (candidate <= 127) return true;

    /* Now do MR - compute Montgomery constants and hash */
    uint64_t n_inv = montgomery_inverse(candidate);
    uint64_t r_sq = montgomery_r_squared(candidate);

    if (!mr_witness_montgomery_cached(candidate, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(candidate, fj64_bases[fj64_hash(candidate)], n_inv, r_sq);
}

static inline uint64_t find_solution_defer_hash(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_defer_hash(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Optimization: Inline Montgomery inverse                                    */
/* ========================================================================== */

static inline uint64_t montgomery_inverse_inline(uint64_t n) {
    /* Newton's method - compute -n^(-1) mod 2^64 */
    uint64_t x = n;  /* Initial approximation (1 bit accurate for odd n) */
    x *= 2 - n * x;  /* 4 bits */
    x *= 2 - n * x;  /* 8 bits */
    x *= 2 - n * x;  /* 16 bits */
    x *= 2 - n * x;  /* 32 bits */
    x *= 2 - n * x;  /* 64 bits */
    return -x;
}

/* For small n, we can use fewer iterations */
static inline uint64_t montgomery_inverse_32bit(uint64_t n) {
    /* For n < 2^32, 4 iterations suffice (32-bit result extended) */
    uint32_t x32 = (uint32_t)n;
    x32 *= 2 - (uint32_t)n * x32;  /* 4 bits */
    x32 *= 2 - (uint32_t)n * x32;  /* 8 bits */
    x32 *= 2 - (uint32_t)n * x32;  /* 16 bits */
    x32 *= 2 - (uint32_t)n * x32;  /* 32 bits */
    /* Extend to 64-bit */
    uint64_t x = x32;
    x *= 2 - n * x;  /* 64 bits */
    return -x;
}

static inline bool is_prime_fj64_fast_inv(uint64_t n) {
    uint64_t n_inv;
    if (n <= UINT32_MAX) {
        n_inv = montgomery_inverse_32bit(n);
    } else {
        n_inv = montgomery_inverse_inline(n);
    }
    uint64_t r_sq = montgomery_r_squared(n);

    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

static inline bool is_candidate_prime_fast_inv(uint64_t candidate) {
    int td = trial_division_current(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast_inv(candidate);
}

static inline uint64_t find_solution_fast_inv(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_fast_inv(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Optimization: Combined best ideas                                          */
/* ========================================================================== */

static inline bool is_candidate_prime_combined(uint64_t candidate) {
    /* Inline first 3 primes (catch ~50%) */
    if (candidate % 3 == 0) return candidate == 3;
    if (candidate % 5 == 0) return candidate == 5;
    if (candidate % 7 == 0) return candidate == 7;

    /* Check rest via loop */
    for (int i = 3; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return candidate == TRIAL_PRIMES[i];
        }
    }

    if (candidate <= 127) return true;

    /* MR with skip-2 optimization */
    uint64_t n_inv = montgomery_inverse(candidate);
    uint64_t r_sq = montgomery_r_squared(candidate);
    uint32_t hash_witness = fj64_bases[fj64_hash(candidate)];

    if (hash_witness == 2) {
        return mr_witness_montgomery_cached(candidate, 2, n_inv, r_sq);
    }

    if (!mr_witness_montgomery_cached(candidate, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(candidate, hash_witness, n_inv, r_sq);
}

static inline uint64_t find_solution_combined(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_combined(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Benchmarks                                                                 */
/* ========================================================================== */

void benchmark(uint64_t n_start, const char* label,
               uint64_t (*find_fn)(uint64_t, uint64_t, uint64_t*)) {
    clock_t start = clock();
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        sum += find_fn(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("%-25s: %.0f n/sec (%.3f sec)\n", label, rate, elapsed);
}

/* Count how often hash witness == 2 */
void count_witness2(uint64_t n_start) {
    uint64_t count2 = 0;
    uint64_t total = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < 1000000ULL; i++) {
        uint64_t candidate = (N - a_max * a_max) >> 1;

        /* Only count candidates that pass trial division */
        bool passes = true;
        for (int j = 0; j < NUM_TRIAL_PRIMES && passes; j++) {
            if (candidate % TRIAL_PRIMES[j] == 0 && candidate != TRIAL_PRIMES[j]) {
                passes = false;
            }
        }

        if (passes && candidate > 127) {
            total++;
            if (fj64_bases[fj64_hash(candidate)] == 2) {
                count2++;
            }
        }

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    printf("Hash witness == 2: %llu / %llu (%.2f%%)\n",
           (unsigned long long)count2, (unsigned long long)total,
           100.0 * count2 / total);
}

int main(void) {
    printf("Optimization Ideas Benchmark Part 3 (%llu iterations)\n\n", ITERATIONS);

    printf("=== Hash witness == 2 frequency ===\n");
    count_witness2(1000000000000ULL);
    count_witness2(2000000000000000000ULL);
    printf("\n");

    printf("=== At n = 10^12 ===\n");
    benchmark(1000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(1000000000000ULL, "Skip base-2 when w=2", find_solution_skip2);
    benchmark(1000000000000ULL, "Defer hash after TD", find_solution_defer_hash);
    benchmark(1000000000000ULL, "Fast 32-bit inverse", find_solution_fast_inv);
    benchmark(1000000000000ULL, "Combined", find_solution_combined);

    printf("\n=== At n = 10^15 ===\n");
    benchmark(1000000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(1000000000000000ULL, "Skip base-2 when w=2", find_solution_skip2);
    benchmark(1000000000000000ULL, "Defer hash after TD", find_solution_defer_hash);
    benchmark(1000000000000000ULL, "Fast 32-bit inverse", find_solution_fast_inv);
    benchmark(1000000000000000ULL, "Combined", find_solution_combined);

    printf("\n=== At n = 2e18 ===\n");
    benchmark(2000000000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(2000000000000000000ULL, "Skip base-2 when w=2", find_solution_skip2);
    benchmark(2000000000000000000ULL, "Defer hash after TD", find_solution_defer_hash);
    benchmark(2000000000000000000ULL, "Fast 32-bit inverse", find_solution_fast_inv);
    benchmark(2000000000000000000ULL, "Combined", find_solution_combined);

    return 0;
}
