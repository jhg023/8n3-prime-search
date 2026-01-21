/*
 * Benchmark trial division with different prime counts
 * Goal: Find optimal number of trial primes for large n values
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "arith.h"
#include "arith_montgomery.h"
#include "fj64_table.h"

/* All 168 primes up to 1000 */
static const uint32_t ALL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,           /* 0-9: 10 primes */
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,        /* 10-19 */
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,  /* 20-29: 30 primes */
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, /* 30-39 */
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233, /* 40-49: 50 primes */
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, /* 50-59 */
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, /* 60-69: 70 primes */
    359, 367, 373, 379, 383, 389, 397, 401, 409, 419, /* 70-79 */
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, /* 80-89: 90 primes */
    479, 487, 491, 499, 503, 509, 521, 523, 541, 547, /* 90-99 */
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, /* 100-109: 110 primes */
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, /* 110-119: 120 primes */
    673, 677, 683, 691, 701, 709, 719, 727, 733, 739, /* 120-129 */
    743, 751, 757, 761, 769, 773, 787, 797, 809, 811, /* 130-139: 140 primes */
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, /* 140-149 */
    881, 883, 887, 907, 911, 919, 929, 937, 941, 947, /* 150-159: 160 primes */
    953, 967, 971, 977, 983, 991, 997, 1009           /* 160-167: 168 primes */
};

/* FJ64 hash */
static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

/* FJ64 with Montgomery */
static inline bool is_prime_fj64_mont(uint64_t n) {
    if (!mr_witness_montgomery(n, 2))
        return false;
    return mr_witness_montgomery(n, fj64_bases[fj64_hash(n)]);
}

/* Configurable trial division + primality test */
static inline bool is_candidate_prime_with_trials(uint64_t candidate, int num_trials) {
    /* Trial division */
    for (int i = 0; i < num_trials; i++) {
        if (candidate % ALL_PRIMES[i] == 0) {
            return (candidate == ALL_PRIMES[i]);
        }
    }

    /* Small primes that passed trial */
    if (candidate <= (uint64_t)ALL_PRIMES[num_trials - 1] * ALL_PRIMES[num_trials - 1]) {
        return true;
    }

    /* Miller-Rabin */
    return is_prime_fj64_mont(candidate);
}

/* Find solution with configurable trial count */
static inline uint64_t find_solution_trials(uint64_t n, int num_trials) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;
        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;
            if (is_candidate_prime_with_trials(candidate, num_trials)) {
                return a;
            }
        }
        if (a < 3) break;
        a -= 2;
    }
    return 0;
}

int main(void) {
    setbuf(stdout, NULL);

    printf("Trial Division Tuning Analysis\n");
    printf("==============================\n\n");

    /* Test configurations */
    int trial_counts[] = {10, 20, 30, 50, 70, 90, 120, 150, 168};
    int num_configs = sizeof(trial_counts) / sizeof(trial_counts[0]);

    /* Test at different scales */
    struct {
        uint64_t n_start;
        uint64_t count;
        const char* label;
    } scales[] = {
        {1000000000000ULL,       500000, "10^12"},
        {1000000000000000ULL,    100000, "10^15"},
        {100000000000000000ULL,   50000, "10^17"},
        {2000000000000000000ULL,  20000, "2e18"}
    };
    int num_scales = sizeof(scales) / sizeof(scales[0]);

    /* Results storage */
    double results[9][4];  /* [config][scale] = n/sec */

    for (int s = 0; s < num_scales; s++) {
        printf("=== Scale: %s (n_start = %llu) ===\n\n",
               scales[s].label, (unsigned long long)scales[s].n_start);

        printf("%-12s  %-12s  %-15s  %-10s\n",
               "Trials", "Max Prime", "Rate (n/sec)", "Relative");
        printf("--------------------------------------------------------\n");

        double baseline = 0;

        for (int c = 0; c < num_configs; c++) {
            int num_trials = trial_counts[c];
            uint64_t n_start = scales[s].n_start;
            uint64_t count = scales[s].count;

            /* Warmup */
            for (uint64_t i = 0; i < 1000; i++) {
                volatile uint64_t a = find_solution_trials(n_start + i, num_trials);
                (void)a;
            }

            /* Timed run */
            clock_t t0 = clock();
            for (uint64_t i = 0; i < count; i++) {
                volatile uint64_t a = find_solution_trials(n_start + i, num_trials);
                (void)a;
            }
            clock_t t1 = clock();

            double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
            double rate = count / elapsed;
            results[c][s] = rate;

            if (c == 0) baseline = rate;

            printf("%-12d  %-12d  %15.0f  %9.2fx\n",
                   num_trials, ALL_PRIMES[num_trials - 1], rate, rate / baseline);
        }

        /* Find optimal */
        int best_config = 0;
        double best_rate = results[0][s];
        for (int c = 1; c < num_configs; c++) {
            if (results[c][s] > best_rate) {
                best_rate = results[c][s];
                best_config = c;
            }
        }
        printf("\nOptimal for %s: %d trials (primes up to %d)\n\n",
               scales[s].label, trial_counts[best_config],
               ALL_PRIMES[trial_counts[best_config] - 1]);
    }

    /* Summary table */
    printf("\n=== Summary: Rate (n/sec) by Configuration ===\n\n");
    printf("%-8s", "Trials");
    for (int s = 0; s < num_scales; s++) {
        printf("  %12s", scales[s].label);
    }
    printf("\n");
    printf("--------");
    for (int s = 0; s < num_scales; s++) {
        printf("  ------------");
    }
    printf("\n");

    for (int c = 0; c < num_configs; c++) {
        printf("%-8d", trial_counts[c]);
        for (int s = 0; s < num_scales; s++) {
            printf("  %12.0f", results[c][s]);
        }
        printf("\n");
    }

    /* Find overall recommendation */
    printf("\n=== Recommendation ===\n\n");

    /* Weight larger scales more heavily since that's the target */
    double weighted_scores[9] = {0};
    double weights[] = {0.1, 0.2, 0.3, 0.4};  /* Weight larger scales more */

    for (int c = 0; c < num_configs; c++) {
        for (int s = 0; s < num_scales; s++) {
            /* Normalize to baseline for that scale */
            weighted_scores[c] += weights[s] * (results[c][s] / results[0][s]);
        }
    }

    int best_overall = 0;
    for (int c = 1; c < num_configs; c++) {
        if (weighted_scores[c] > weighted_scores[best_overall]) {
            best_overall = c;
        }
    }

    printf("For large n optimization (weighted toward 10^17-2e18):\n");
    printf("  Recommended: %d trial primes (up to %d)\n",
           trial_counts[best_overall], ALL_PRIMES[trial_counts[best_overall] - 1]);
    printf("  Current:     120 trial primes (up to 661)\n");

    return 0;
}
