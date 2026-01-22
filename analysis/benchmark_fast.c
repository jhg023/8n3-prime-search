/*
 * Benchmark: Compare baseline vs micro-optimized implementation
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <inttypes.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/prime.h"
#include "../include/solve.h"
#include "../include/solve_fast.h"

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

typedef uint64_t (*solver_fn)(uint64_t n, uint64_t* p_out);

static void benchmark(const char* name, solver_fn fn, uint64_t start_n, uint64_t count) {
    double t0 = get_time();

    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t p;
        uint64_t a = fn(n, &p);

        /* Verify */
        if (a > 0) {
            uint64_t N = 8 * n + 3;
            if (a * a + 2 * p != N) {
                printf("ERROR in %s: n=%" PRIu64 "\n", name, n);
                exit(1);
            }
        }
    }

    double elapsed = get_time() - t0;
    printf("  %-30s %10.0f n/sec  (%.3fs)\n", name, count / elapsed, elapsed);
}

int main(int argc, char** argv) {
    printf("Fast Implementation Benchmark\n");
    printf("=============================\n");

    uint64_t count = 200000;
    if (argc > 1) count = strtoull(argv[1], NULL, 10);

    /* Verify correctness */
    printf("\nVerifying correctness (n = 1 to 10000)...\n");
    for (uint64_t n = 1; n <= 10000; n++) {
        uint64_t N = 8 * n + 3;

        uint64_t p1, p2, p3;
        uint64_t a1 = find_solution(n, &p1);
        uint64_t a2 = find_solution_fast(n, &p2);
        uint64_t a3 = find_solution_fast_v2(n, &p3);

        if (a1 == 0 || a2 == 0 || a3 == 0) {
            printf("ERROR: No solution for n=%" PRIu64 "\n", n);
            return 1;
        }

        if (a1*a1 + 2*p1 != N || a2*a2 + 2*p2 != N || a3*a3 + 2*p3 != N) {
            printf("ERROR: Invalid solution for n=%" PRIu64 "\n", n);
            return 1;
        }
    }
    printf("All implementations correct.\n");

    /* Benchmark at multiple scales */
    uint64_t scales[] = {
        1000000000ULL,
        1000000000000ULL,
        1000000000000000ULL,
        1000000000000000000ULL,
        2000000000000000000ULL,
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        printf("\nBenchmark at n = %.2e, count = %" PRIu64 "\n", (double)scales[i], count);
        printf("%-32s %15s\n", "Implementation", "Throughput");
        printf("--------------------------------------------------------\n");

        benchmark("Baseline (solve.h)", find_solution, scales[i], count);
        benchmark("Fast (solve_fast.h)", find_solution_fast, scales[i], count);
        benchmark("Fast v2 (no incremental)", find_solution_fast_v2, scales[i], count);
    }

    return 0;
}
