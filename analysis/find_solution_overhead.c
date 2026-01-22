/*
 * Measure the overhead of find_solution() vs a version with pre-computed a_max
 *
 * Tests if passing (N, a_max) directly saves measurable time
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define SOLVE_TRACK_STATS
#include "../include/solve.h"

#define ITERATIONS 10000000

/* Version with pre-computed N and a_max */
static inline uint64_t find_solution_from_N(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    /* Initialize candidate = (N - a_max²) / 2 and delta = 2*(a_max - 1) */
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    /* Iterate in reverse: largest a first means smallest candidate p first */
    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime(candidate)) {
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

/* Benchmark original find_solution */
void benchmark_original(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        uint64_t a = find_solution(n_start + i, &p);
        sum += a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("Original find_solution(): %.2f n/sec (sum=%llu)\n",
           rate, (unsigned long long)sum);
}

/* Benchmark with pre-computed N and a_max */
void benchmark_precomputed(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    /* Initialize */
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        uint64_t a = find_solution_from_N(N, a_max, &p);
        sum += a;

        /* Update for next n */
        N += 8;
        uint64_t next_sq = (a_max + 2) * (a_max + 2);
        if (next_sq <= N) {
            a_max += 2;
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("Precomputed (N, a_max):   %.2f n/sec (sum=%llu)\n",
           rate, (unsigned long long)sum);
}

int main(int argc, char **argv) {
    uint64_t n_start = 1000000000000ULL;
    if (argc > 1) {
        n_start = strtoull(argv[1], NULL, 10);
    }

    printf("Testing find_solution optimization at n_start = %llu\n", (unsigned long long)n_start);
    printf("Iterations: %d\n\n", ITERATIONS);

    printf("=== n = 10^12 ===\n");
    benchmark_original(n_start);
    benchmark_precomputed(n_start);

    printf("\n=== n = 10^15 ===\n");
    n_start = 1000000000000000ULL;
    benchmark_original(n_start);
    benchmark_precomputed(n_start);

    printf("\n=== n = 2e18 ===\n");
    n_start = 2000000000000000000ULL;
    benchmark_original(n_start);
    benchmark_precomputed(n_start);

    return 0;
}
