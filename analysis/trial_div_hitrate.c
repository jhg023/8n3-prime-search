/*
 * Analyze trial division hit rates to optimize the prime list.
 * Goal: Find the minimum set of primes that gives maximum filter rate.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <inttypes.h>

#include "../include/arith.h"

/* Extended prime list for analysis */
static const uint32_t ALL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233
};
#define NUM_ALL_PRIMES (sizeof(ALL_PRIMES) / sizeof(ALL_PRIMES[0]))

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char** argv) {
    uint64_t start_n = 1000000000000000ULL;  /* 10^15 */
    uint64_t count = 100000;

    if (argc > 1) start_n = strtoull(argv[1], NULL, 10);
    if (argc > 2) count = strtoull(argv[2], NULL, 10);

    printf("Trial Division Hit Rate Analysis\n");
    printf("================================\n\n");
    printf("Start n: %.2e, Count: %" PRIu64 "\n\n", (double)start_n, count);

    /* Track hits per prime */
    uint64_t hits[NUM_ALL_PRIMES] = {0};
    uint64_t total_candidates = 0;
    uint64_t total_filtered = 0;

    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t N = 8 * n + 3;
        uint64_t a_max = isqrt64(N);
        if ((a_max & 1) == 0) a_max--;

        /* Iterate through candidates */
        for (uint64_t a = a_max; a >= 1; a -= 2) {
            uint64_t a_sq = a * a;
            if (a_sq > N - 4) continue;

            uint64_t candidate = (N - a_sq) >> 1;
            total_candidates++;

            /* Check each prime */
            bool filtered = false;
            for (size_t i = 0; i < NUM_ALL_PRIMES; i++) {
                if (candidate % ALL_PRIMES[i] == 0) {
                    hits[i]++;
                    if (candidate != ALL_PRIMES[i]) {
                        filtered = true;
                    }
                    break;  /* Count first hit only */
                }
            }
            if (filtered) total_filtered++;

            /* Early exit when we find a prime (simulating real search) */
            /* Actually, for this analysis we want to see all candidates */
        }
    }

    printf("Total candidates examined: %" PRIu64 "\n", total_candidates);
    printf("Filtered by trial division: %" PRIu64 " (%.2f%%)\n\n",
           total_filtered, 100.0 * total_filtered / total_candidates);

    printf("Hit rate per prime (first hit):\n");
    printf("%-8s %12s %10s %10s\n", "Prime", "Hits", "Hit Rate", "Cumulative");
    printf("--------------------------------------------\n");

    uint64_t cumulative = 0;
    for (size_t i = 0; i < NUM_ALL_PRIMES; i++) {
        cumulative += hits[i];
        double hit_rate = 100.0 * hits[i] / total_candidates;
        double cum_rate = 100.0 * cumulative / total_candidates;
        printf("%-8" PRIu32 " %12" PRIu64 " %9.2f%% %9.2f%%\n",
               ALL_PRIMES[i], hits[i], hit_rate, cum_rate);
    }

    printf("\n");
    printf("Candidates surviving trial division: %" PRIu64 " (%.2f%%)\n",
           total_candidates - total_filtered,
           100.0 * (total_candidates - total_filtered) / total_candidates);

    /* Analyze optimal cutoffs */
    printf("\nOptimal prime count analysis:\n");
    printf("%-10s %12s %10s\n", "Primes", "Filter Rate", "Marginal");
    printf("----------------------------------\n");

    cumulative = 0;
    double prev_rate = 0;
    for (size_t i = 0; i < NUM_ALL_PRIMES; i++) {
        cumulative += hits[i];
        double rate = 100.0 * cumulative / total_candidates;
        double marginal = rate - prev_rate;
        if ((i+1) % 5 == 0 || i < 10) {
            printf("%-10zu %11.2f%% %9.2f%%\n", i+1, rate, marginal);
        }
        prev_rate = rate;
    }

    return 0;
}
