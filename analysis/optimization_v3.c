/*
 * Optimization v3: Precomputed residue tables
 *
 * Key insight: For N = 8n + 3 and candidate p = (N - a²)/2,
 * whether p is divisible by small prime q depends only on:
 *   - n mod q
 *   - a mod 2q
 *
 * We can precompute a bitmap of "bad" (n mod M, a mod 2M) pairs
 * for M = LCM(3, 5, 7) = 105. This eliminates per-N computation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/prime.h"
#include "../include/solve.h"

/* ========================================================================== */
/* Timing utilities                                                           */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Precomputed tables                                                         */
/* ========================================================================== */

/*
 * For primes 3, 5, 7: M = 105, 2M = 210
 * We need to check all (n mod 105, a mod 210) pairs.
 * A pair is "bad" if the candidate is divisible by 3, 5, or 7.
 *
 * Table size: 105 * 105 = 11,025 entries (since odd a only, we have 105 odd residues mod 210)
 * Actually, a mod 210 where a is odd: 1, 3, 5, ..., 209 = 105 values
 *
 * We use a bitmap: 105 rows, each with 105 bits packed into uint64_t[2]
 */

#define TABLE_N_MOD 105
#define TABLE_A_MOD 210  /* But we only track odd a, so 105 values */
#define NUM_ODD_RESIDUES 105

/* Bitmap: skip_table[n_mod][word] where each bit represents an odd a residue */
static uint64_t skip_table[TABLE_N_MOD][2];  /* 105 bits need 2 uint64_t */
static bool table_initialized = false;

static void init_skip_table(void) {
    if (table_initialized) return;

    memset(skip_table, 0, sizeof(skip_table));

    for (int n_mod = 0; n_mod < TABLE_N_MOD; n_mod++) {
        /* N = 8n + 3, so N mod 210 depends on n mod 105 and position within 210 period */
        /* But actually we need N mod (2*q) for each q to check divisibility */

        /* For each odd a_mod in [0, 210) */
        for (int a_idx = 0; a_idx < NUM_ODD_RESIDUES; a_idx++) {
            int a_mod = 2 * a_idx + 1;  /* Odd residues: 1, 3, 5, ..., 209 */

            /* Check if candidate is divisible by 3, 5, or 7 */
            /* candidate = (N - a²) / 2 = (8n + 3 - a²) / 2 */
            /* We need to check if (8n + 3 - a²) / 2 ≡ 0 (mod q) */
            /* i.e., 8n + 3 - a² ≡ 0 (mod 2q) */

            bool should_skip = false;

            /* Check mod 3 */
            {
                int n3 = n_mod % 3;
                int a3 = a_mod % 6;  /* a mod 2*3 */
                int a_sq_mod6 = (a3 * a3) % 6;
                int N_mod6 = (8 * n3 + 3) % 6;
                if ((N_mod6 - a_sq_mod6 + 6) % 6 == 0) {
                    /* candidate divisible by 3, but check if it equals 3 */
                    /* For small candidates this could be 3, but at scale it's rare */
                    should_skip = true;
                }
            }

            /* Check mod 5 */
            if (!should_skip) {
                int n5 = n_mod % 5;
                int a5 = a_mod % 10;
                int a_sq_mod10 = (a5 * a5) % 10;
                int N_mod10 = (8 * n5 + 3) % 10;
                if ((N_mod10 - a_sq_mod10 + 10) % 10 == 0) {
                    should_skip = true;
                }
            }

            /* Check mod 7 */
            if (!should_skip) {
                int n7 = n_mod % 7;
                int a7 = a_mod % 14;
                int a_sq_mod14 = (a7 * a7) % 14;
                int N_mod14 = (8 * n7 + 3) % 14;
                if ((N_mod14 - a_sq_mod14 + 14) % 14 == 0) {
                    should_skip = true;
                }
            }

            if (should_skip) {
                /* Set bit a_idx in skip_table[n_mod] */
                int word = a_idx / 64;
                int bit = a_idx % 64;
                skip_table[n_mod][word] |= (1ULL << bit);
            }
        }
    }

    table_initialized = true;

    /* Print statistics */
    int total_skip = 0;
    for (int n_mod = 0; n_mod < TABLE_N_MOD; n_mod++) {
        total_skip += __builtin_popcountll(skip_table[n_mod][0]) +
                      __builtin_popcountll(skip_table[n_mod][1]);
    }
    /* Silent - uncomment for debugging:
    printf("Skip table: %d / %d pairs marked (%.1f%%)\n",
           total_skip, TABLE_N_MOD * NUM_ODD_RESIDUES,
           100.0 * total_skip / (TABLE_N_MOD * NUM_ODD_RESIDUES));
    */
}

static inline bool should_skip_lookup(int n_mod, int a_idx) {
    int word = a_idx / 64;
    int bit = a_idx % 64;
    return (skip_table[n_mod][word] >> bit) & 1;
}

/* ========================================================================== */
/* Strategy: Precomputed table lookup                                         */
/* ========================================================================== */

static uint64_t find_solution_table(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    int n_mod = n % TABLE_N_MOD;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;

        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;

            /* For small candidates (3, 5, 7), always check - table may incorrectly skip */
            if (candidate <= 7) {
                if (candidate == 3 || candidate == 5 || candidate == 7) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            } else {
                /* Check skip table */
                int a_mod = a % TABLE_A_MOD;
                int a_idx = a_mod / 2;

                if (!should_skip_lookup(n_mod, a_idx)) {
                    if (is_candidate_prime(candidate)) {
                        if (p_out) *p_out = candidate;
                        return a;
                    }
                }
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Strategy: Table + reduced trial division                                    */
/* ========================================================================== */

/* Trial primes excluding 3, 5, 7 (handled by table) */
static const uint32_t REDUCED_PRIMES[] = {
    11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
    47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    97, 101, 103, 107, 109, 113, 127
};
#define NUM_REDUCED_PRIMES (sizeof(REDUCED_PRIMES) / sizeof(REDUCED_PRIMES[0]))

static inline bool is_prime_table_reduced(uint64_t candidate) {
    /* Trial division (3, 5, 7 already filtered by table) */
    for (size_t i = 0; i < NUM_REDUCED_PRIMES; i++) {
        if (candidate % REDUCED_PRIMES[i] == 0) {
            return candidate == REDUCED_PRIMES[i];
        }
    }
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_table_reduced(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    int n_mod = n % TABLE_N_MOD;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;

        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;

            /* Handle small primes that table might skip incorrectly */
            if (candidate <= 7) {
                if (candidate == 3 || candidate == 5 || candidate == 7) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            } else {
                int a_mod = a % TABLE_A_MOD;
                int a_idx = a_mod / 2;

                if (!should_skip_lookup(n_mod, a_idx)) {
                    if (is_prime_table_reduced(candidate)) {
                        if (p_out) *p_out = candidate;
                        return a;
                    }
                }
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Baseline (reference)                                                       */
/* ========================================================================== */

static uint64_t find_solution_baseline(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;

        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;

            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Benchmarking                                                               */
/* ========================================================================== */

typedef uint64_t (*solver_fn)(uint64_t n, uint64_t* p_out);

typedef struct {
    const char* name;
    solver_fn solve;
} Strategy;

static Strategy strategies[] = {
    {"Baseline", find_solution_baseline},
    {"Precomputed table", find_solution_table},
    {"Table + reduced trial", find_solution_table_reduced},
};

#define NUM_STRATEGIES (sizeof(strategies) / sizeof(strategies[0]))

static void benchmark_strategy(Strategy* s, uint64_t start_n, uint64_t count) {
    double t0 = get_time();
    uint64_t solutions = 0;

    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t p;
        uint64_t a = s->solve(n, &p);
        if (a > 0) {
            solutions++;
            uint64_t N = 8 * n + 3;
            if (a * a + 2 * p != N) {
                printf("ERROR: Verification failed for n=%" PRIu64 "\n", n);
                exit(1);
            }
        }
    }

    double elapsed = get_time() - t0;
    double rate = count / elapsed;

    printf("  %-28s %10.0f n/sec  (%.3fs)\n", s->name, rate, elapsed);
}

static void run_benchmarks(uint64_t start_n, uint64_t count) {
    printf("\nBenchmark at n = %.2e, count = %" PRIu64 "\n", (double)start_n, count);
    printf("%-30s %15s\n", "Strategy", "Throughput");
    printf("------------------------------------------------------\n");

    for (size_t i = 0; i < NUM_STRATEGIES; i++) {
        benchmark_strategy(&strategies[i], start_n, count);
    }
}

int main(int argc, char** argv) {
    printf("Optimization v3: Precomputed Residue Tables\n");
    printf("============================================\n");

    init_skip_table();

    uint64_t count = 100000;
    if (argc > 1) count = strtoull(argv[1], NULL, 10);

    /* Verify correctness */
    printf("\nVerifying correctness (n = 1 to 10000)...\n");
    for (uint64_t n = 1; n <= 10000; n++) {
        uint64_t N = 8 * n + 3;

        for (size_t i = 0; i < NUM_STRATEGIES; i++) {
            uint64_t p;
            uint64_t a = strategies[i].solve(n, &p);

            if (a == 0) {
                printf("ERROR: %s found no solution for n=%" PRIu64 "\n",
                       strategies[i].name, n);
                return 1;
            }

            if (a * a + 2 * p != N) {
                printf("ERROR: %s invalid solution for n=%" PRIu64 ": a=%" PRIu64 " p=%" PRIu64 "\n",
                       strategies[i].name, n, a, p);
                return 1;
            }
        }
    }
    printf("All strategies produce valid solutions.\n");

    /* Benchmark */
    uint64_t scales[] = {
        1000000000ULL,
        1000000000000ULL,
        1000000000000000ULL,
        1000000000000000000ULL,
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        run_benchmarks(scales[i], count);
    }

    return 0;
}
