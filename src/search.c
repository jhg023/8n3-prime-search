/*
 * Counterexample Search: 8n + 3 = a^2 + 2p
 *
 * Searches for counterexamples to the conjecture that every integer
 * of the form 8n + 3 can be expressed as a^2 + 2p, where a is a
 * positive odd integer and p is prime.
 *
 * Uses FJ64_262K primality test for optimal performance:
 * - Only 2 Miller-Rabin tests (vs 7 in standard deterministic test)
 * - 512KB precomputed hash table for witness selection
 * - 100% deterministic for all 64-bit integers
 *
 * Reference: Forisek & Jancina (2015), "Fast Primality Testing for
 * Integers That Fit into a Machine Word"
 *
 * Compile: make release
 * Usage:   ./search [n_start] [n_end]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

/* Enable statistics tracking in solve.h */
#define SOLVE_TRACK_STATS

/* Include shared headers */
#include "fmt.h"
#include "solve.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Progress reporting interval in seconds */
#define PROGRESS_SECONDS 5.0

/* Default search range */
#define DEFAULT_N_START 1000000000000ULL    /* 10^12 */
#define DEFAULT_N_END   1000010000000ULL    /* 10^12 + 10^7 */

/* ========================================================================== */
/* Search Loop                                                                */
/* ========================================================================== */

/**
 * Run the search over a range of n values
 *
 * Uses incremental N and a_max tracking for better performance:
 * - N = 8n + 3 increases by 8 each step
 * - a_max = isqrt(N) changes very rarely (every ~a_max/2 iterations)
 */
void run_search(uint64_t n_start, uint64_t n_end,
                uint64_t *out_counterexamples) {
    clock_t start_time = clock();

    uint64_t counterexamples_found = 0;
    double last_progress_time = 0.0;

    /* Reset statistics */
    solve_reset_stats();

    /* Initialize N and a_max for incremental tracking */
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t n = n_start; n < n_end; n++) {
        uint64_t p;
        uint64_t a = find_solution_from_N(N, a_max, &p);

        if (a == 0) {
            /* Counterexample found! */
            counterexamples_found++;
            printf("COUNTEREXAMPLE: n = %s (N = %s)\n",
                   fmt_num(n), fmt_num(N));
            fflush(stdout);
        }

        /* Update N and a_max for next iteration */
        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) {
            a_max = next_a;
        }

        /* Progress reporting - check every 262144 iterations (~6 checks/sec at 1.5M n/sec) */
        if ((n & 0x3FFFF) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= PROGRESS_SECONDS) {
                uint64_t n_done, total_checks, bits32;
                solve_get_stats(&n_done, &total_checks, &bits32);
                double avg_checks = solve_get_avg_checks();
                double pct_32bit = (total_checks > 0) ? 100.0 * bits32 / total_checks : 0.0;

                double rate = (n - n_start + 1) / elapsed;
                double pct = 100.0 * (n - n_start) / (n_end - n_start);

                printf("n = %s (%.1f%%), rate = %s n/sec, avg_checks = %.2f, 32-bit: %.1f%%\n",
                       fmt_num(n), pct, fmt_num((uint64_t)rate), avg_checks, pct_32bit);
                fflush(stdout);

                last_progress_time = elapsed;
            }
        }
    }

    *out_counterexamples = counterexamples_found;
}

/* ========================================================================== */
/* Verification                                                               */
/* ========================================================================== */

/**
 * Verify the algorithm against known solutions
 */
bool verify_known_solutions(void) {
    printf("Verifying known solutions...\n");

    struct { uint64_t n; uint64_t a; uint64_t p; } known[] = {
        {1, 1, 5},    /* 11 = 1 + 10 */
        {2, 3, 5},    /* 19 = 9 + 10 */
        {3, 1, 13},   /* 27 = 1 + 26 */
        {4, 5, 5}     /* 35 = 25 + 10 (also 1 + 34 = 1 + 2*17) */
    };

    bool all_pass = true;

    for (int i = 0; i < 4; i++) {
        uint64_t n = known[i].n;
        uint64_t expected_a = known[i].a;
        uint64_t expected_p = known[i].p;

        uint64_t N = 8 * n + 3;
        uint64_t lhs = N;
        uint64_t rhs = expected_a * expected_a + 2 * expected_p;

        bool equation_valid = (lhs == rhs);
        bool p_is_prime = is_prime_64(expected_p);

        uint64_t found_p;
        uint64_t found_a = find_solution(n, &found_p);

        printf("  n=%llu: N=%llu, given (%llu,%llu), found (%llu,%llu) ... ",
               (unsigned long long)n, (unsigned long long)N,
               (unsigned long long)expected_a, (unsigned long long)expected_p,
               (unsigned long long)found_a, (unsigned long long)found_p);

        if (equation_valid && p_is_prime && found_a > 0) {
            printf("PASS\n");
        } else {
            printf("FAIL\n");
            all_pass = false;
        }
    }

    return all_pass;
}

/* ========================================================================== */
/* Argument Parsing                                                           */
/* ========================================================================== */

/**
 * Parse a number that may be in scientific notation (e.g., 1e12, 2.5e9)
 */
uint64_t parse_number(const char* str) {
    char* endptr;
    double val = strtod(str, &endptr);
    if (*endptr != '\0') {
        return strtoull(str, NULL, 10);
    }
    return (uint64_t)val;
}

void print_usage(const char* program) {
    printf("Usage: %s [n_start] [n_end]\n", program);
    printf("\n");
    printf("Search for counterexamples to: 8n + 3 = a^2 + 2p\n");
    printf("\n");
    printf("For each n in [n_start, n_end), attempts to find odd a and prime p\n");
    printf("such that 8n + 3 = a^2 + 2p. Reports any n for which no solution exists.\n");
    printf("\n");
    printf("Arguments:\n");
    printf("  n_start   Starting value of n (inclusive), default: 1e12\n");
    printf("  n_end     Ending value of n (exclusive), default: 1e12 + 1e7\n");
    printf("\n");
    printf("Numbers can be in scientific notation (e.g., 1e9, 2.5e6)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s                    Search [10^12, 10^12 + 10^7)\n", program);
    printf("  %s 1 1e6              Search [1, 10^6)\n", program);
    printf("  %s 1e9 2e9            Search [10^9, 2*10^9)\n", program);
    printf("\n");
    printf("Exit codes:\n");
    printf("  0  Search completed, no counterexamples found\n");
    printf("  1  Error (invalid arguments, verification failure)\n");
    printf("  2  Counterexample found\n");
}

/* ========================================================================== */
/* Main                                                                       */
/* ========================================================================== */

int main(int argc, char** argv) {
    /* Disable stdout buffering for real-time output */
    setbuf(stdout, NULL);

    printf("==================================================================\n");
    printf("     Counterexample Search: 8n + 3 = a^2 + 2p                     \n");
    printf("     Optimized with FJ64_262K primality test                      \n");
    printf("==================================================================\n\n");

    uint64_t n_start = DEFAULT_N_START;
    uint64_t n_end = DEFAULT_N_END;

    /* Handle help flag */
    if (argc == 2 && (strcmp(argv[1], "-h") == 0 ||
                      strcmp(argv[1], "--help") == 0)) {
        print_usage(argv[0]);
        return 0;
    }

    /* Parse arguments */
    if (argc >= 2) {
        n_start = parse_number(argv[1]);
    }
    if (argc >= 3) {
        n_end = parse_number(argv[2]);
    }

    if (n_start >= n_end) {
        fprintf(stderr, "Error: n_start must be less than n_end\n");
        return 1;
    }

    uint64_t total = n_end - n_start;

    printf("Configuration:\n");
    printf("  Range: n in [%s, %s)\n", fmt_num(n_start), fmt_num(n_end));
    printf("  Count: %s values\n", fmt_num(total));
    printf("  Primality test: FJ64_262K (2 Miller-Rabin tests)\n");
    printf("\n");

    /* Verify algorithm correctness */
    if (!verify_known_solutions()) {
        fprintf(stderr, "\nERROR: Verification failed!\n");
        return 1;
    }
    printf("\n");

    /* Run search */
    printf("Starting search...\n\n");
    clock_t global_start = clock();

    uint64_t total_counterexamples = 0;

    run_search(n_start, n_end, &total_counterexamples);

    clock_t global_end = clock();
    double global_elapsed = (double)(global_end - global_start) / CLOCKS_PER_SEC;

    /* Get final statistics */
    uint64_t stat_n, stat_checks, stat_32bit;
    solve_get_stats(&stat_n, &stat_checks, &stat_32bit);
    double avg_checks = solve_get_avg_checks();
    double pct_32bit = (stat_checks > 0) ? 100.0 * stat_32bit / stat_checks : 0.0;

    /* Print results */
    printf("\n");
    printf("==================================================================\n");
    printf("RESULTS\n");
    printf("==================================================================\n\n");

    printf("Total time:           %.1f seconds\n", global_elapsed);
    printf("Total throughput:     %s n/sec\n",
           fmt_num((uint64_t)(total / global_elapsed)));
    printf("Counterexamples:      %s\n", fmt_num(total_counterexamples));
    printf("Avg checks per n:     %.2f\n", avg_checks);
    printf("Total a's checked:    %s\n", fmt_num(stat_checks));
    printf("32-bit candidates:    %s (%.2f%%)\n",
           fmt_num(stat_32bit), pct_32bit);

    return (total_counterexamples > 0) ? 2 : 0;
}
