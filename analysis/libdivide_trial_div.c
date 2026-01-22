/*
 * Benchmark: libdivide for trial division optimization
 *
 * Trial division currently uses candidate % prime for 30 small primes.
 * libdivide can replace division with multiply+shift operations.
 *
 * Key insight: We need MODULO (remainder), not quotient.
 * libdivide gives us quotient q = n / d
 * We compute remainder as: r = n - q * d
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

/* Trial primes (same as prime.h) */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_TRIAL_PRIMES 30

#define ITERATIONS 10000000ULL

/* Pre-generated libdivide structures for each trial prime */
static struct libdivide_u64_t trial_dividers[NUM_TRIAL_PRIMES];
static struct libdivide_u64_branchfree_t trial_dividers_bf[NUM_TRIAL_PRIMES];

/* Initialize libdivide structures */
void init_libdivide(void) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        trial_dividers[i] = libdivide_u64_gen(TRIAL_PRIMES[i]);
        trial_dividers_bf[i] = libdivide_u64_branchfree_gen(TRIAL_PRIMES[i]);
    }
}

/* ========================================================================== */
/* Original trial division (using % operator)                                 */
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
/* libdivide trial division (branchfull)                                      */
/* ========================================================================== */

static inline int trial_division_libdivide(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        uint64_t q = libdivide_u64_do(candidate, &trial_dividers[i]);
        uint64_t r = candidate - q * TRIAL_PRIMES[i];
        if (r == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* ========================================================================== */
/* libdivide trial division (branchfree)                                      */
/* ========================================================================== */

static inline int trial_division_libdivide_bf(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        uint64_t q = libdivide_u64_branchfree_do(candidate, &trial_dividers_bf[i]);
        uint64_t r = candidate - q * TRIAL_PRIMES[i];
        if (r == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* ========================================================================== */
/* Benchmarks                                                                 */
/* ========================================================================== */

/* Generate test candidates (mix of primes and composites) */
void generate_candidates(uint64_t *candidates, uint64_t n_start, uint64_t count) {
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < count; i++) {
        candidates[i] = (N - a_max * a_max) >> 1;
        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }
}

void benchmark_original(uint64_t *candidates, uint64_t count) {
    clock_t start = clock();
    volatile int sum = 0;

    for (uint64_t i = 0; i < count; i++) {
        sum += trial_division_original(candidates[i]);
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double ns_per_call = (elapsed * 1e9) / count;

    printf("Original (%%):         %.2f ns/call, %.3f sec (sum=%d)\n",
           ns_per_call, elapsed, sum);
}

void benchmark_libdivide(uint64_t *candidates, uint64_t count) {
    clock_t start = clock();
    volatile int sum = 0;

    for (uint64_t i = 0; i < count; i++) {
        sum += trial_division_libdivide(candidates[i]);
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double ns_per_call = (elapsed * 1e9) / count;

    printf("libdivide (branchfull): %.2f ns/call, %.3f sec (sum=%d)\n",
           ns_per_call, elapsed, sum);
}

void benchmark_libdivide_bf(uint64_t *candidates, uint64_t count) {
    clock_t start = clock();
    volatile int sum = 0;

    for (uint64_t i = 0; i < count; i++) {
        sum += trial_division_libdivide_bf(candidates[i]);
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double ns_per_call = (elapsed * 1e9) / count;

    printf("libdivide (branchfree): %.2f ns/call, %.3f sec (sum=%d)\n",
           ns_per_call, elapsed, sum);
}

/* Verify correctness */
void verify_correctness(uint64_t *candidates, uint64_t count) {
    int mismatches = 0;
    for (uint64_t i = 0; i < count && i < 100000; i++) {
        int orig = trial_division_original(candidates[i]);
        int ld = trial_division_libdivide(candidates[i]);
        int ld_bf = trial_division_libdivide_bf(candidates[i]);

        if (orig != ld || orig != ld_bf) {
            mismatches++;
            if (mismatches <= 5) {
                printf("MISMATCH at candidate=%llu: orig=%d, ld=%d, ld_bf=%d\n",
                       (unsigned long long)candidates[i], orig, ld, ld_bf);
            }
        }
    }
    printf("Verification: %d mismatches\n\n", mismatches);
}

int main(void) {
    printf("Benchmarking libdivide for trial division (%llu iterations)\n\n", ITERATIONS);

    /* Initialize libdivide */
    init_libdivide();

    /* Allocate and generate candidates */
    uint64_t *candidates = malloc(ITERATIONS * sizeof(uint64_t));
    if (!candidates) {
        fprintf(stderr, "Failed to allocate memory\n");
        return 1;
    }

    printf("=== At n = 10^12 (small candidates) ===\n");
    generate_candidates(candidates, 1000000000000ULL, ITERATIONS);
    verify_correctness(candidates, ITERATIONS);
    benchmark_original(candidates, ITERATIONS);
    benchmark_libdivide(candidates, ITERATIONS);
    benchmark_libdivide_bf(candidates, ITERATIONS);

    printf("\n=== At n = 10^15 (medium candidates) ===\n");
    generate_candidates(candidates, 1000000000000000ULL, ITERATIONS);
    verify_correctness(candidates, ITERATIONS);
    benchmark_original(candidates, ITERATIONS);
    benchmark_libdivide(candidates, ITERATIONS);
    benchmark_libdivide_bf(candidates, ITERATIONS);

    printf("\n=== At n = 2e18 (large candidates) ===\n");
    generate_candidates(candidates, 2000000000000000000ULL, ITERATIONS);
    verify_correctness(candidates, ITERATIONS);
    benchmark_original(candidates, ITERATIONS);
    benchmark_libdivide(candidates, ITERATIONS);
    benchmark_libdivide_bf(candidates, ITERATIONS);

    free(candidates);
    return 0;
}
