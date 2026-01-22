/*
 * Test progress checking with simple counter vs bitmask
 */

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "../include/solve.h"

#define ITERATIONS 10000000ULL
#define PROGRESS_INTERVAL 16384

/* Current: bitmask on n */
void benchmark_bitmask(uint64_t n_start) {
    clock_t start_time = clock();
    double last_progress_time = 0.0;
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t n = n_start; n < n_start + ITERATIONS; n++) {
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;

        if ((n & 0x3FFF) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= 5.0) {
                last_progress_time = elapsed;
            }
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("Bitmask (n & 0x3FFF):  %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

/* Alternative: simple counter */
void benchmark_counter(uint64_t n_start) {
    clock_t start_time = clock();
    double last_progress_time = 0.0;
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint32_t progress_counter = PROGRESS_INTERVAL;

    for (uint64_t n = n_start; n < n_start + ITERATIONS; n++) {
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;

        if (--progress_counter == 0) {
            progress_counter = PROGRESS_INTERVAL;
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= 5.0) {
                last_progress_time = elapsed;
            }
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("Counter (--counter):   %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

/* Alternative: iteration index check */
void benchmark_index(uint64_t n_start) {
    clock_t start_time = clock();
    double last_progress_time = 0.0;
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;

        if ((i & 0x3FFF) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= 5.0) {
                last_progress_time = elapsed;
            }
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("Index (i & 0x3FFF):    %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

/* No progress at all */
void benchmark_none(uint64_t n_start) {
    clock_t start_time = clock();
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("No progress check:     %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

int main(void) {
    printf("=== At n=2e18 ===\n");
    uint64_t n_start = 2000000000000000000ULL;
    benchmark_none(n_start);
    benchmark_bitmask(n_start);
    benchmark_index(n_start);
    benchmark_counter(n_start);

    printf("\n=== At n=10^12 ===\n");
    n_start = 1000000000000ULL;
    benchmark_none(n_start);
    benchmark_bitmask(n_start);
    benchmark_index(n_start);
    benchmark_counter(n_start);

    return 0;
}
