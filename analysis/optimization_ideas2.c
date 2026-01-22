/*
 * Optimization Ideas Benchmark - Part 2
 *
 * Test combinations and more aggressive optimizations:
 * 1. Combined: unrolled + no check
 * 2. Fully unrolled trial division (all 30 primes)
 * 3. GCD-based filtering (check gcd with primorial)
 * 4. Early termination tweaks in Montgomery MR
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
/* Optimization: Fully unrolled trial division (all 30)                       */
/* ========================================================================== */

static inline int trial_division_full_unroll(uint64_t c) {
    /* Fully unrolled - no loop overhead at all */
    if (c % 3 == 0) return (c == 3) ? 1 : 0;
    if (c % 5 == 0) return (c == 5) ? 1 : 0;
    if (c % 7 == 0) return (c == 7) ? 1 : 0;
    if (c % 11 == 0) return (c == 11) ? 1 : 0;
    if (c % 13 == 0) return (c == 13) ? 1 : 0;
    if (c % 17 == 0) return (c == 17) ? 1 : 0;
    if (c % 19 == 0) return (c == 19) ? 1 : 0;
    if (c % 23 == 0) return (c == 23) ? 1 : 0;
    if (c % 29 == 0) return (c == 29) ? 1 : 0;
    if (c % 31 == 0) return (c == 31) ? 1 : 0;
    if (c % 37 == 0) return (c == 37) ? 1 : 0;
    if (c % 41 == 0) return (c == 41) ? 1 : 0;
    if (c % 43 == 0) return (c == 43) ? 1 : 0;
    if (c % 47 == 0) return (c == 47) ? 1 : 0;
    if (c % 53 == 0) return (c == 53) ? 1 : 0;
    if (c % 59 == 0) return (c == 59) ? 1 : 0;
    if (c % 61 == 0) return (c == 61) ? 1 : 0;
    if (c % 67 == 0) return (c == 67) ? 1 : 0;
    if (c % 71 == 0) return (c == 71) ? 1 : 0;
    if (c % 73 == 0) return (c == 73) ? 1 : 0;
    if (c % 79 == 0) return (c == 79) ? 1 : 0;
    if (c % 83 == 0) return (c == 83) ? 1 : 0;
    if (c % 89 == 0) return (c == 89) ? 1 : 0;
    if (c % 97 == 0) return (c == 97) ? 1 : 0;
    if (c % 101 == 0) return (c == 101) ? 1 : 0;
    if (c % 103 == 0) return (c == 103) ? 1 : 0;
    if (c % 107 == 0) return (c == 107) ? 1 : 0;
    if (c % 109 == 0) return (c == 109) ? 1 : 0;
    if (c % 113 == 0) return (c == 113) ? 1 : 0;
    if (c % 127 == 0) return (c == 127) ? 1 : 0;
    return 2;
}

static inline bool is_candidate_prime_full_unroll(uint64_t candidate) {
    int td = trial_division_full_unroll(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_full_unroll(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_full_unroll(candidate)) {
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
/* Optimization: Split small primes (3,5,7) vs rest                           */
/* ========================================================================== */

/* These catch ~50% of composites */
static inline int trial_division_split(uint64_t c) {
    /* Check 3, 5, 7 first (catch ~50% of composites) */
    if (c % 3 == 0) return (c == 3) ? 1 : 0;
    if (c % 5 == 0) return (c == 5) ? 1 : 0;
    if (c % 7 == 0) return (c == 7) ? 1 : 0;

    /* Then check 11-127 via loop */
    for (int i = 3; i < NUM_TRIAL_PRIMES; i++) {
        if (c % TRIAL_PRIMES[i] == 0) {
            return (c == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_split(uint64_t candidate) {
    int td = trial_division_split(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_split(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_split(candidate)) {
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
/* Optimization: GCD with primorial                                           */
/* ========================================================================== */

/* Product of first 8 primes = 3*5*7*11*13*17*19*23 = 223092870 */
#define PRIMORIAL_8 223092870ULL

static inline uint64_t gcd64(uint64_t a, uint64_t b) {
    while (b) {
        uint64_t t = b;
        b = a % b;
        a = t;
    }
    return a;
}

static inline int trial_division_gcd(uint64_t c) {
    /* First check GCD with primorial of first 8 primes */
    uint64_t g = gcd64(c, PRIMORIAL_8);
    if (g > 1) {
        /* c is divisible by one of 3,5,7,11,13,17,19,23 */
        /* Check which one */
        if (c % 3 == 0) return (c == 3) ? 1 : 0;
        if (c % 5 == 0) return (c == 5) ? 1 : 0;
        if (c % 7 == 0) return (c == 7) ? 1 : 0;
        if (c % 11 == 0) return (c == 11) ? 1 : 0;
        if (c % 13 == 0) return (c == 13) ? 1 : 0;
        if (c % 17 == 0) return (c == 17) ? 1 : 0;
        if (c % 19 == 0) return (c == 19) ? 1 : 0;
        if (c % 23 == 0) return (c == 23) ? 1 : 0;
    }

    /* Check remaining primes 29-127 */
    for (int i = 8; i < NUM_TRIAL_PRIMES; i++) {
        if (c % TRIAL_PRIMES[i] == 0) {
            return (c == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_gcd(uint64_t candidate) {
    int td = trial_division_gcd(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_gcd(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_gcd(candidate)) {
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
/* Optimization: Simplified trial (just first 8 primes)                       */
/* ========================================================================== */

static inline int trial_division_8only(uint64_t c) {
    /* Only check first 8 primes - ~75% filter rate */
    if (c % 3 == 0) return (c == 3) ? 1 : 0;
    if (c % 5 == 0) return (c == 5) ? 1 : 0;
    if (c % 7 == 0) return (c == 7) ? 1 : 0;
    if (c % 11 == 0) return (c == 11) ? 1 : 0;
    if (c % 13 == 0) return (c == 13) ? 1 : 0;
    if (c % 17 == 0) return (c == 17) ? 1 : 0;
    if (c % 19 == 0) return (c == 19) ? 1 : 0;
    if (c % 23 == 0) return (c == 23) ? 1 : 0;
    return 2;
}

static inline bool is_candidate_prime_8only(uint64_t candidate) {
    int td = trial_division_8only(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 23) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_8only(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_8only(candidate)) {
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
    printf("Optimization Ideas Benchmark Part 2 (%llu iterations)\n\n", ITERATIONS);

    printf("=== At n = 10^12 ===\n");
    benchmark(1000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(1000000000000ULL, "Full unroll (30)", find_solution_full_unroll);
    benchmark(1000000000000ULL, "Split 3,5,7 + loop", find_solution_split);
    benchmark(1000000000000ULL, "GCD primorial", find_solution_gcd);
    benchmark(1000000000000ULL, "Only 8 primes", find_solution_8only);

    printf("\n=== At n = 10^15 ===\n");
    benchmark(1000000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(1000000000000000ULL, "Full unroll (30)", find_solution_full_unroll);
    benchmark(1000000000000000ULL, "Split 3,5,7 + loop", find_solution_split);
    benchmark(1000000000000000ULL, "GCD primorial", find_solution_gcd);
    benchmark(1000000000000000ULL, "Only 8 primes", find_solution_8only);

    printf("\n=== At n = 2e18 ===\n");
    benchmark(2000000000000000000ULL, "Current (baseline)", find_solution_current);
    benchmark(2000000000000000000ULL, "Full unroll (30)", find_solution_full_unroll);
    benchmark(2000000000000000000ULL, "Split 3,5,7 + loop", find_solution_split);
    benchmark(2000000000000000000ULL, "GCD primorial", find_solution_gcd);
    benchmark(2000000000000000000ULL, "Only 8 primes", find_solution_8only);

    return 0;
}
