/*
 * Measure isqrt64 overhead in the main search loop
 *
 * Compare:
 * 1. Current: isqrt64(N) called per-n
 * 2. Optimized: a_max maintained incrementally
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include "../include/arith.h"

#define ITERATIONS 10000000

/* Measure isqrt64 cost per call */
void benchmark_isqrt64(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t N = 8 * (n_start + i) + 3;
        sum += isqrt64(N);
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double ns_per_call = (elapsed * 1e9) / ITERATIONS;

    printf("isqrt64 per-n:     %.2f ns/call (sum=%llu to prevent optimize-out)\n",
           ns_per_call, (unsigned long long)sum);
}

/* Measure incremental a_max tracking */
void benchmark_incremental(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        /* Use a_max */
        sum += a_max;

        /* Increment N */
        N += 8;

        /* Occasionally adjust a_max (when (a_max+2)^2 <= N) */
        uint64_t next_sq = (a_max + 2) * (a_max + 2);
        if (next_sq <= N) {
            a_max += 2;
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double ns_per_call = (elapsed * 1e9) / ITERATIONS;

    printf("Incremental a_max: %.2f ns/iter (sum=%llu)\n",
           ns_per_call, (unsigned long long)sum);
}

/* Verify correctness of incremental approach */
void verify_incremental(uint64_t n_start, uint64_t count) {
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max_inc = isqrt64(N);
    if ((a_max_inc & 1) == 0) a_max_inc--;

    int mismatches = 0;

    for (uint64_t i = 0; i < count; i++) {
        /* Compute expected a_max */
        uint64_t N_cur = 8 * (n_start + i) + 3;
        uint64_t a_max_expected = isqrt64(N_cur);
        if ((a_max_expected & 1) == 0) a_max_expected--;

        if (a_max_inc != a_max_expected) {
            mismatches++;
            if (mismatches <= 5) {
                printf("MISMATCH at n=%llu: expected %llu, got %llu\n",
                       (unsigned long long)(n_start + i),
                       (unsigned long long)a_max_expected,
                       (unsigned long long)a_max_inc);
            }
        }

        /* Increment for next iteration */
        N += 8;
        uint64_t next_sq = (a_max_inc + 2) * (a_max_inc + 2);
        if (next_sq <= N) {
            a_max_inc += 2;
        }
    }

    printf("Verification: %d mismatches out of %llu\n", mismatches, (unsigned long long)count);
}

/* Count how often a_max changes */
void count_amax_changes(uint64_t n_start, uint64_t count) {
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t changes = 0;

    for (uint64_t i = 0; i < count; i++) {
        N += 8;
        uint64_t next_sq = (a_max + 2) * (a_max + 2);
        if (next_sq <= N) {
            a_max += 2;
            changes++;
        }
    }

    printf("a_max changes: %llu times over %llu iterations (every %.1f iterations on average)\n",
           (unsigned long long)changes, (unsigned long long)count,
           (double)count / (changes > 0 ? changes : 1));
}

int main(int argc, char **argv) {
    uint64_t n_start = 1000000000000ULL; /* 10^12 */
    if (argc > 1) {
        n_start = strtoull(argv[1], NULL, 10);
    }

    printf("Testing isqrt64 optimization at n_start = %llu\n\n", (unsigned long long)n_start);

    /* First verify correctness */
    printf("=== Verification ===\n");
    verify_incremental(n_start, 100000);
    printf("\n");

    /* Count a_max changes */
    printf("=== a_max change frequency ===\n");
    count_amax_changes(n_start, ITERATIONS);
    printf("\n");

    /* Benchmark */
    printf("=== Benchmark (%d iterations) ===\n", ITERATIONS);
    benchmark_isqrt64(n_start);
    benchmark_incremental(n_start);

    /* Also test at 2e18 scale */
    printf("\n=== At n = 2e18 scale ===\n");
    n_start = 2000000000000000000ULL;
    count_amax_changes(n_start, ITERATIONS);
    benchmark_isqrt64(n_start);
    benchmark_incremental(n_start);

    return 0;
}
