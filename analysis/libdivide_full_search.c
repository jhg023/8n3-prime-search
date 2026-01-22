/*
 * Benchmark: libdivide in full search context
 *
 * Test if libdivide helps when integrated into the actual find_solution loop
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#include "libdivide.h"
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

/* Pre-generated libdivide structures */
static struct libdivide_u64_branchfree_t trial_dividers[NUM_TRIAL_PRIMES];

void init_libdivide(void) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        trial_dividers[i] = libdivide_u64_branchfree_gen(TRIAL_PRIMES[i]);
    }
}

/* ========================================================================== */
/* Original trial division                                                    */
/* ========================================================================== */

static inline int trial_division_original(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* ========================================================================== */
/* libdivide trial division                                                   */
/* ========================================================================== */

static inline int trial_division_libdivide(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        uint64_t q = libdivide_u64_branchfree_do(candidate, &trial_dividers[i]);
        uint64_t r = candidate - q * TRIAL_PRIMES[i];
        if (r == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* ========================================================================== */
/* FJ64 primality (shared)                                                    */
/* ========================================================================== */

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline bool is_prime_fj64_fast_local(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

/* ========================================================================== */
/* Full find_solution implementations                                         */
/* ========================================================================== */

static inline bool is_candidate_prime_original(uint64_t candidate) {
    int td = trial_division_original(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast_local(candidate);
}

static inline bool is_candidate_prime_libdivide(uint64_t candidate) {
    int td = trial_division_libdivide(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast_local(candidate);
}

static inline uint64_t find_solution_original(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_original(candidate)) {
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

static inline uint64_t find_solution_libdivide(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_libdivide(candidate)) {
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

    printf("%-20s: %.0f n/sec, %.3f sec\n", label, rate, elapsed);
}

int main(void) {
    printf("Benchmarking libdivide in full search context (%llu iterations)\n\n", ITERATIONS);

    init_libdivide();

    printf("=== n = 10^12 ===\n");
    benchmark(1000000000000ULL, "Original", find_solution_original);
    benchmark(1000000000000ULL, "libdivide", find_solution_libdivide);

    printf("\n=== n = 10^15 ===\n");
    benchmark(1000000000000000ULL, "Original", find_solution_original);
    benchmark(1000000000000000ULL, "libdivide", find_solution_libdivide);

    printf("\n=== n = 2e18 ===\n");
    benchmark(2000000000000000000ULL, "Original", find_solution_original);
    benchmark(2000000000000000000ULL, "libdivide", find_solution_libdivide);

    return 0;
}
