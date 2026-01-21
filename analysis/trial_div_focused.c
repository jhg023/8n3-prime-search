/*
 * Focused comparison: 30 vs 50 vs 70 vs 120 trial primes
 * More iterations for statistical significance
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "arith.h"
#include "arith_montgomery.h"
#include "fj64_table.h"

/* Primes for trial division */
static const uint32_t PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,           /* 10 */
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,        /* 20 */
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,  /* 30 */
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, /* 40 */
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233, /* 50 */
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, /* 60 */
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, /* 70 */
    359, 367, 373, 379, 383, 389, 397, 401, 409, 419, /* 80 */
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, /* 90 */
    479, 487, 491, 499, 503, 509, 521, 523, 541, 547, /* 100 */
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607, /* 110 */
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661  /* 120 */
};

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline bool is_prime_fj64_mont(uint64_t n) {
    if (!mr_witness_montgomery(n, 2))
        return false;
    return mr_witness_montgomery(n, fj64_bases[fj64_hash(n)]);
}

/* Macro to generate specialized functions for each trial count */
#define MAKE_SOLVER(NUM_TRIALS) \
static inline uint64_t find_solution_##NUM_TRIALS(uint64_t n) { \
    uint64_t N = 8 * n + 3; \
    uint64_t a_max = isqrt64(N); \
    if ((a_max & 1) == 0) a_max--; \
    uint64_t a = a_max; \
    while (1) { \
        uint64_t a_sq = a * a; \
        if (a_sq <= N - 4) { \
            uint64_t candidate = (N - a_sq) >> 1; \
            /* Trial division */ \
            bool filtered = false; \
            for (int i = 0; i < NUM_TRIALS; i++) { \
                if (candidate % PRIMES[i] == 0) { \
                    if (candidate == PRIMES[i]) return a; \
                    filtered = true; \
                    break; \
                } \
            } \
            if (!filtered) { \
                if (candidate <= 127) return a; \
                if (is_prime_fj64_mont(candidate)) return a; \
            } \
        } \
        if (a < 3) break; \
        a -= 2; \
    } \
    return 0; \
}

MAKE_SOLVER(30)
MAKE_SOLVER(50)
MAKE_SOLVER(70)
MAKE_SOLVER(120)

typedef uint64_t (*solver_fn)(uint64_t);

int main(void) {
    setbuf(stdout, NULL);

    printf("Focused Trial Division Comparison\n");
    printf("==================================\n\n");

    struct {
        int num_trials;
        solver_fn solver;
    } configs[] = {
        {30, find_solution_30},
        {50, find_solution_50},
        {70, find_solution_70},
        {120, find_solution_120}
    };
    int num_configs = 4;

    struct {
        uint64_t n_start;
        uint64_t count;
        const char* label;
    } scales[] = {
        {1000000000000ULL,       1000000, "10^12"},
        {1000000000000000ULL,     200000, "10^15"},
        {100000000000000000ULL,   100000, "10^17"},
        {2000000000000000000ULL,   50000, "2e18"}
    };
    int num_scales = 4;

    printf("Running with larger iteration counts for accuracy...\n\n");

    for (int s = 0; s < num_scales; s++) {
        printf("=== %s ===\n", scales[s].label);
        printf("%-8s  %12s  %10s\n", "Trials", "Rate (n/sec)", "vs 120");
        printf("--------------------------------\n");

        double rate_120 = 0;

        for (int c = 0; c < num_configs; c++) {
            uint64_t n_start = scales[s].n_start;
            uint64_t count = scales[s].count;
            solver_fn solver = configs[c].solver;

            /* Warmup */
            for (uint64_t i = 0; i < 10000; i++) {
                volatile uint64_t a = solver(n_start + i);
                (void)a;
            }

            /* Timed run */
            clock_t t0 = clock();
            for (uint64_t i = 0; i < count; i++) {
                volatile uint64_t a = solver(n_start + i);
                (void)a;
            }
            clock_t t1 = clock();

            double elapsed = (double)(t1 - t0) / CLOCKS_PER_SEC;
            double rate = count / elapsed;

            if (configs[c].num_trials == 120) rate_120 = rate;

            const char* marker = "";
            if (c == num_configs - 1) {
                /* Find best */
                double best_rate = rate;
                int best_idx = c;
                for (int j = 0; j < c; j++) {
                    /* Re-run to get rate (hacky but works) */
                }
            }

            printf("%-8d  %12.0f  %+9.1f%%\n",
                   configs[c].num_trials, rate,
                   (rate_120 > 0) ? 100.0 * (rate - rate_120) / rate_120 : 0.0);
        }
        printf("\n");
    }

    return 0;
}
