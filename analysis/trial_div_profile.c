/*
 * Profile trial division: which primes are doing the filtering?
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>

#include "arith.h"

/* All 120 primes */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,           /* 0-9 */
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,        /* 10-19 */
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,  /* 20-29 */
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, /* 30-39 */
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233, /* 40-49 */
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, /* 50-59 */
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, /* 60-69 */
    359, 367, 373, 379, 383, 389, 397, 401, 409, 419, /* 70-79 */
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, /* 80-89 */
    479, 487, 491, 499, 503, 509, 521, 523, 541, 547, /* 90-99 */
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, /* 100-109 */
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661  /* 110-119 */
};

static uint64_t hits[120] = {0};

int main(void) {
    printf("Trial Division Hit Rate Analysis\n");
    printf("=================================\n\n");

    /* Test at 2e18 scale */
    uint64_t n_start = 2000000000000000000ULL;
    uint64_t count = 10000;

    printf("Testing at n ~ 2e18 with %llu n values\n\n", (unsigned long long)count);

    uint64_t total_candidates = 0;
    uint64_t total_filtered = 0;
    uint64_t passed_all = 0;

    for (uint64_t i = 0; i < count; i++) {
        uint64_t n = n_start + i;
        uint64_t N = 8 * n + 3;
        uint64_t a_max = isqrt64(N);
        if ((a_max & 1) == 0) a_max--;

        for (uint64_t a = a_max; a >= 1; a -= 2) {
            uint64_t a_sq = a * a;
            if (a_sq > N - 4) continue;
            uint64_t candidate = (N - a_sq) >> 1;
            total_candidates++;

            bool filtered = false;
            for (int j = 0; j < 120; j++) {
                if (candidate % TRIAL_PRIMES[j] == 0) {
                    hits[j]++;
                    filtered = true;
                    total_filtered++;
                    break;  /* Only count first hit */
                }
            }
            if (!filtered) passed_all++;

            /* Early exit if we found a probable prime (simplified) */
            if (!filtered && candidate > 127) break;
        }
    }

    printf("Summary:\n");
    printf("  Total candidates:     %llu\n", (unsigned long long)total_candidates);
    printf("  Filtered by trial:    %llu (%.1f%%)\n",
           (unsigned long long)total_filtered, 100.0 * total_filtered / total_candidates);
    printf("  Passed all trials:    %llu (%.1f%%)\n",
           (unsigned long long)passed_all, 100.0 * passed_all / total_candidates);

    printf("\nHit rate by prime (first hit only):\n");
    printf("  Prime   Hits       Cumulative%%\n");
    printf("  ------  ---------  -----------\n");

    uint64_t cumulative = 0;
    for (int i = 0; i < 120; i++) {
        cumulative += hits[i];
        if (hits[i] > 0 || i < 30 || i == 119) {
            printf("  %3d     %9llu  %5.1f%%\n",
                   TRIAL_PRIMES[i], (unsigned long long)hits[i],
                   100.0 * cumulative / total_candidates);
        }
    }

    /* Find cutoff points */
    printf("\nCumulative filter rate at various cutoffs:\n");
    int cutoffs[] = {10, 20, 30, 50, 70, 100, 120};
    for (int c = 0; c < 7; c++) {
        uint64_t sum = 0;
        for (int i = 0; i < cutoffs[c]; i++) sum += hits[i];
        printf("  First %3d primes (up to %3d): %.2f%% filtered\n",
               cutoffs[c], TRIAL_PRIMES[cutoffs[c]-1], 100.0 * sum / total_candidates);
    }

    return 0;
}
