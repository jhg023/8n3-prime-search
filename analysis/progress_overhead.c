/*
 * Measure overhead of different progress checking strategies
 */

#include <stdio.h>
#include <stdint.h>
#include <time.h>

#define ITERATIONS 100000000ULL

/* Measure clock() call overhead */
void benchmark_clock_calls(void) {
    clock_t start = clock();
    volatile clock_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        if ((i & 0x3FFF) == 0) {
            sum += clock();
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    uint64_t num_calls = ITERATIONS / 16384;
    double ns_per_call = (elapsed * 1e9) / num_calls;

    printf("clock() every 16384 iters: %.2f ns/call, %llu calls, %.3f sec total\n",
           ns_per_call, (unsigned long long)num_calls, elapsed);
}

/* Measure overhead of just the bitmask check */
void benchmark_no_clock(void) {
    clock_t start = clock();
    volatile uint64_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        if ((i & 0x3FFF) == 0) {
            sum += i;  /* Do something trivial instead of clock() */
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Bitmask only (no clock): %.3f sec (sum=%llu)\n", elapsed, (unsigned long long)sum);
}

/* Measure overhead with larger check interval */
void benchmark_larger_interval(void) {
    clock_t start = clock();
    volatile clock_t sum = 0;

    /* Check every ~1M iterations (0xFFFFF = 1048575) */
    for (uint64_t i = 0; i < ITERATIONS; i++) {
        if ((i & 0xFFFFF) == 0) {
            sum += clock();
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    uint64_t num_calls = ITERATIONS / (0xFFFFF + 1);

    printf("clock() every 1M iters:   %.3f sec, %llu calls\n",
           elapsed, (unsigned long long)num_calls);
}

/* Counter-based approach: decrement counter, check clock when zero */
void benchmark_counter_approach(void) {
    clock_t start = clock();
    volatile clock_t sum = 0;

    uint64_t check_counter = 1000000;  /* Check every 1M iterations */

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        if (--check_counter == 0) {
            sum += clock();
            check_counter = 1000000;
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Counter decrement:        %.3f sec\n", elapsed);
}

int main(void) {
    printf("Testing progress check overhead (%llu iterations)\n\n", ITERATIONS);

    benchmark_no_clock();
    benchmark_clock_calls();
    benchmark_larger_interval();
    benchmark_counter_approach();

    return 0;
}
