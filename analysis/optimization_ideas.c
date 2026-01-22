/*
 * Optimization Ideas Benchmark
 *
 * Test various micro-optimizations identified in code review:
 * 1. Trial division loop unrolling (4x)
 * 2. Removing candidate >= 2 check for large n
 * 3. Simplified FJ64 hash
 * 4. Combined optimizations
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

static inline uint32_t fj64_hash_current(uint64_t x) {
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
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash_current(n)], n_inv, r_sq);
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
/* Optimization 1: Unrolled trial division (4x)                               */
/* ========================================================================== */

static inline int trial_division_unrolled(uint64_t candidate) {
    int i = 0;

    /* Unroll 4x for first 28 primes */
    for (; i < 28; i += 4) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
        if (candidate % TRIAL_PRIMES[i+1] == 0) {
            return (candidate == TRIAL_PRIMES[i+1]) ? 1 : 0;
        }
        if (candidate % TRIAL_PRIMES[i+2] == 0) {
            return (candidate == TRIAL_PRIMES[i+2]) ? 1 : 0;
        }
        if (candidate % TRIAL_PRIMES[i+3] == 0) {
            return (candidate == TRIAL_PRIMES[i+3]) ? 1 : 0;
        }
    }

    /* Handle remaining 2 primes */
    if (candidate % TRIAL_PRIMES[28] == 0) {
        return (candidate == TRIAL_PRIMES[28]) ? 1 : 0;
    }
    if (candidate % TRIAL_PRIMES[29] == 0) {
        return (candidate == TRIAL_PRIMES[29]) ? 1 : 0;
    }

    return 2;
}

static inline bool is_candidate_prime_unrolled(uint64_t candidate) {
    int td = trial_division_unrolled(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_unrolled(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_unrolled(candidate)) {
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
/* Optimization 2: Remove candidate >= 2 check (for large n)                  */
/* ========================================================================== */

/* At n >= 10^9, the minimum candidate p = (N - a_max^2)/2 is always > 2 */
static inline uint64_t find_solution_no_check(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        /* Skip the candidate >= 2 check - always true for large n */
        if (is_candidate_prime_current(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Optimization 3: Simplified FJ64 hash (fewer ops)                           */
/* ========================================================================== */

/* Try single-round hash (likely loses correctness, but test speed) */
static inline uint32_t fj64_hash_simple(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

/* Two-round hash (still needs verification) */
static inline uint32_t fj64_hash_two(uint64_t x) {
    x ^= x >> 32;
    x *= 0x45d9f3b3335b369ULL;
    x ^= x >> 32;
    return x & 262143;
}

/* We can't actually use these without verifying the hash table works with them */

/* ========================================================================== */
/* Optimization 4: Inline the most common path                                */
/* ========================================================================== */

/* Inline trial division for first 8 primes (most likely to hit) */
static inline bool is_candidate_prime_inline8(uint64_t candidate) {
    /* Inline first 8 primes (catch ~75% of composites) */
    if (candidate % 3 == 0) return candidate == 3;
    if (candidate % 5 == 0) return candidate == 5;
    if (candidate % 7 == 0) return candidate == 7;
    if (candidate % 11 == 0) return candidate == 11;
    if (candidate % 13 == 0) return candidate == 13;
    if (candidate % 17 == 0) return candidate == 17;
    if (candidate % 19 == 0) return candidate == 19;
    if (candidate % 23 == 0) return candidate == 23;

    /* Rest via loop */
    for (int i = 8; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]);
        }
    }

    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_inline8(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_inline8(candidate)) {
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

int main(void) {
    printf("Optimization Ideas Benchmark (%llu iterations)\n\n", ITERATIONS);

    printf("=== At n = 10^12 ===\n");
    benchmark(1000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(1000000000000ULL, "Unrolled trial div (4x)", find_solution_unrolled);
    benchmark(1000000000000ULL, "No candidate>=2 check", find_solution_no_check);
    benchmark(1000000000000ULL, "Inline first 8 primes", find_solution_inline8);

    printf("\n=== At n = 10^15 ===\n");
    benchmark(1000000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(1000000000000000ULL, "Unrolled trial div (4x)", find_solution_unrolled);
    benchmark(1000000000000000ULL, "No candidate>=2 check", find_solution_no_check);
    benchmark(1000000000000000ULL, "Inline first 8 primes", find_solution_inline8);

    printf("\n=== At n = 2e18 ===\n");
    benchmark(2000000000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(2000000000000000000ULL, "Unrolled trial div (4x)", find_solution_unrolled);
    benchmark(2000000000000000000ULL, "No candidate>=2 check", find_solution_no_check);
    benchmark(2000000000000000000ULL, "Inline first 8 primes", find_solution_inline8);

    return 0;
}
