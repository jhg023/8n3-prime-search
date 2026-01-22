/*
 * Test different progress check intervals
 */

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#include "../include/solve.h"

#define ITERATIONS 10000000ULL

void benchmark_interval(uint64_t n_start, uint32_t mask, const char* label) {
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

        if ((i & mask) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= 5.0) {
                last_progress_time = elapsed;
            }
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start_time) / CLOCKS_PER_SEC;
    uint64_t num_checks = ITERATIONS / (mask + 1);
    printf("%-20s: %.3f sec, %.0f n/sec (%llu clock() calls)\n",
           label, elapsed, ITERATIONS / elapsed, (unsigned long long)num_checks);
}

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
    printf("%-20s: %.3f sec, %.0f n/sec (0 clock() calls)\n",
           "No check", elapsed, ITERATIONS / elapsed);
}

int main(void) {
    printf("=== At n=2e18 (%llu iterations) ===\n", ITERATIONS);
    uint64_t n_start = 2000000000000000000ULL;

    benchmark_none(n_start);
    benchmark_interval(n_start, 0x3FFF,   "16K (0x3FFF)");
    benchmark_interval(n_start, 0xFFFF,   "64K (0xFFFF)");
    benchmark_interval(n_start, 0x3FFFF,  "256K (0x3FFFF)");
    benchmark_interval(n_start, 0xFFFFF,  "1M (0xFFFFF)");

    printf("\n=== At n=10^12 ===\n");
    n_start = 1000000000000ULL;

    benchmark_none(n_start);
    benchmark_interval(n_start, 0x3FFF,   "16K (0x3FFF)");
    benchmark_interval(n_start, 0xFFFF,   "64K (0xFFFF)");
    benchmark_interval(n_start, 0x3FFFF,  "256K (0x3FFFF)");
    benchmark_interval(n_start, 0xFFFFF,  "1M (0xFFFFF)");

    return 0;
}
