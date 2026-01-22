/*
 * Measure overhead of progress checking in realistic search context
 */

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "../include/solve.h"

#define ITERATIONS 10000000ULL

/* Current implementation: bitmask + clock() + floating point */
void benchmark_current(uint64_t n_start) {
    clock_t start_time = clock();
    double last_progress_time = 0.0;
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t n = n_start + i;
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;

        /* Current progress check */
        if ((n & 0x3FFF) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= 5.0) {
                last_progress_time = elapsed;
                /* Would print here, but skip for benchmark */
            }
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("Current (16K mask + clock): %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

/* No progress checking at all (baseline) */
void benchmark_no_progress(uint64_t n_start) {
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
    printf("No progress check:          %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

/* Larger interval: check every ~1M iterations */
void benchmark_larger_interval(uint64_t n_start) {
    clock_t start_time = clock();
    double last_progress_time = 0.0;
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t n = n_start + i;
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;

        /* Larger interval: 0x7FFFF = 524287 (~500K) */
        if ((n & 0x7FFFF) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= 5.0) {
                last_progress_time = elapsed;
            }
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("Larger interval (512K):     %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

/* Counter-based: estimate iterations per 5 seconds */
void benchmark_counter_based(uint64_t n_start) {
    clock_t start_time = clock();
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    /* Initial estimate: assume 1M n/sec, so 5M iterations per check */
    uint64_t check_interval = 5000000;
    uint64_t next_check = check_interval;
    double last_progress_time = 0.0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        sum += find_solution_from_N(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;

        /* Counter-based check */
        if (i >= next_check) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;

            if (elapsed - last_progress_time >= 5.0) {
                /* Adjust interval based on actual rate */
                double rate = (i + 1) / elapsed;
                check_interval = (uint64_t)(rate * 5.0);  /* ~5 seconds worth */
                if (check_interval < 100000) check_interval = 100000;
                last_progress_time = elapsed;
            }

            next_check = i + check_interval;
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    printf("Counter-based (adaptive):   %.3f sec, %.0f n/sec\n",
           elapsed, ITERATIONS / elapsed);
}

int main(void) {
    uint64_t n_start = 2000000000000000000ULL;  /* 2e18 - worst case */

    printf("Testing progress check strategies at n=2e18 (%llu iterations)\n\n",
           ITERATIONS);

    benchmark_no_progress(n_start);
    benchmark_current(n_start);
    benchmark_larger_interval(n_start);
    benchmark_counter_based(n_start);

    printf("\n=== At n=10^12 ===\n");
    n_start = 1000000000000ULL;
    benchmark_no_progress(n_start);
    benchmark_current(n_start);
    benchmark_larger_interval(n_start);
    benchmark_counter_based(n_start);

    return 0;
}
