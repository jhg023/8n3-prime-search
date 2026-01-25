/*
 * Batched Counterexample Search: 8n + 3 = a^2 + 2p
 *
 * Uses segmented sieve approach to process batches of n values together,
 * exploiting the arithmetic progression structure of candidate primes.
 *
 * This approach is optimized for exhaustive verification (checking ALL n
 * in a range) rather than per-n early exit.
 *
 * Compile: make search_batched
 * Usage:   ./search_batched [n_start] [n_end] [--batch-size N]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* Include shared headers */
#include "fmt.h"
#include "arith.h"
#include "prime.h"
#include "batch_sieve.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Progress reporting interval in seconds */
#define PROGRESS_SECONDS 5.0

/* Default search range */
#define DEFAULT_N_START 1000000000ULL      /* 10^9 */
#define DEFAULT_N_END   1010000000ULL      /* 10^9 + 10^7 */

/* Default batch size */
#define DEFAULT_BATCH_SIZE 65536

/* ========================================================================== */
/* Time Formatting                                                            */
/* ========================================================================== */

static const char* fmt_time(double seconds) {
    static char buf[64];

    if (seconds < 60) {
        snprintf(buf, sizeof(buf), "%.0fs", seconds);
    } else if (seconds < 3600) {
        int mins = (int)(seconds / 60);
        int secs = (int)(seconds) % 60;
        snprintf(buf, sizeof(buf), "%dm %ds", mins, secs);
    } else if (seconds < 86400) {
        int hours = (int)(seconds / 3600);
        int mins = ((int)(seconds) % 3600) / 60;
        snprintf(buf, sizeof(buf), "%dh %dm", hours, mins);
    } else {
        int days = (int)(seconds / 86400);
        int hours = ((int)(seconds) % 86400) / 3600;
        snprintf(buf, sizeof(buf), "%dd %dh", days, hours);
    }
    return buf;
}

/* ========================================================================== */
/* Argument Parsing                                                           */
/* ========================================================================== */

uint64_t parse_number(const char* str) {
    char* endptr;
    double val = strtod(str, &endptr);
    if (*endptr != '\0') {
        return strtoull(str, NULL, 10);
    }
    return (uint64_t)val;
}

void print_usage(const char* program) {
    printf("Usage: %s [n_start] [n_end] [options]\n", program);
    printf("\n");
    printf("Batched search for counterexamples to: 8n + 3 = a^2 + 2p\n");
    printf("\n");
    printf("Uses segmented sieve to process batches of n values together.\n");
    printf("Optimized for exhaustive verification of large ranges.\n");
    printf("\n");
    printf("Arguments:\n");
    printf("  n_start            Starting value of n (inclusive), default: 1e9\n");
    printf("  n_end              Ending value of n (exclusive), default: 1e9 + 1e7\n");
    printf("  --batch-size N     Batch size (default: 65536)\n");
    printf("\n");
    printf("Numbers can be in scientific notation (e.g., 1e9, 2.5e6)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s                        Search [10^9, 10^9 + 10^7)\n", program);
    printf("  %s 1 1e6                  Search [1, 10^6)\n", program);
    printf("  %s 1e9 2e9 --batch-size 131072  Use larger batches\n", program);
    printf("\n");
    printf("Exit codes:\n");
    printf("  0  Search completed, no counterexamples found\n");
    printf("  1  Error (invalid arguments)\n");
    printf("  2  Counterexample found\n");
}

/* ========================================================================== */
/* Verification                                                               */
/* ========================================================================== */

/**
 * Verify a single solution
 */
bool verify_solution(uint64_t n, uint64_t a, uint64_t p) {
    uint64_t N = 8 * n + 3;
    uint64_t computed = a * a + 2 * p;
    return (N == computed) && (a & 1) && is_prime_64(p);
}

/* ========================================================================== */
/* Main                                                                       */
/* ========================================================================== */

int main(int argc, char** argv) {
    /* Disable stdout buffering for real-time output */
    setbuf(stdout, NULL);

    printf("==================================================================\n");
    printf("     Batched Counterexample Search: 8n + 3 = a^2 + 2p            \n");
    printf("     Segmented sieve approach for exhaustive verification        \n");
    printf("==================================================================\n\n");

    uint64_t n_start = DEFAULT_N_START;
    uint64_t n_end = DEFAULT_N_END;
    uint64_t batch_size = DEFAULT_BATCH_SIZE;

    /* Handle help flag */
    if (argc >= 2 && (strcmp(argv[1], "-h") == 0 ||
                      strcmp(argv[1], "--help") == 0)) {
        print_usage(argv[0]);
        return 0;
    }

    /* Parse arguments */
    int arg_idx = 1;
    int pos_count = 0;
    while (arg_idx < argc) {
        if (strcmp(argv[arg_idx], "--batch-size") == 0 && arg_idx + 1 < argc) {
            batch_size = parse_number(argv[arg_idx + 1]);
            arg_idx += 2;
        } else if (argv[arg_idx][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[arg_idx]);
            print_usage(argv[0]);
            return 1;
        } else {
            /* Positional argument */
            if (pos_count == 0) {
                n_start = parse_number(argv[arg_idx]);
            } else if (pos_count == 1) {
                n_end = parse_number(argv[arg_idx]);
            }
            pos_count++;
            arg_idx++;
        }
    }

    /* Set default end if only start was given */
    if (pos_count == 1) {
        n_end = n_start + 10000000;
    }

    if (n_start >= n_end) {
        fprintf(stderr, "Error: n_start must be less than n_end\n");
        return 1;
    }

    if (batch_size < 1024) {
        fprintf(stderr, "Warning: batch_size too small, using 1024\n");
        batch_size = 1024;
    }

    uint64_t total = n_end - n_start;
    uint64_t num_batches = (total + batch_size - 1) / batch_size;

    printf("Configuration:\n");
    printf("  Range: n in [%s, %s)\n", fmt_num(n_start), fmt_num(n_end));
    printf("  Count: %s values\n", fmt_num(total));
    printf("  Batch size: %s\n", fmt_num(batch_size));
    printf("  Number of batches: %s\n", fmt_num(num_batches));
    printf("\n");

    /* Create batch sieve context */
    BatchSieve *bs = batch_sieve_create(n_start, batch_size);
    if (!bs) {
        fprintf(stderr, "Error: Failed to allocate batch sieve\n");
        return 1;
    }

    /* Statistics */
    uint64_t total_solved = 0;
    uint64_t total_mr_saved = 0;
    uint64_t total_mr_done = 0;
    uint64_t total_counterexamples = 0;

    /* Run search */
    printf("Starting batched search...\n\n");

    double global_start;
#ifdef _OPENMP
    global_start = omp_get_wtime();
#else
    global_start = (double)clock() / CLOCKS_PER_SEC;
#endif

    double last_report_time = 0.0;
    uint64_t batches_processed = 0;

    for (uint64_t batch_start = n_start; batch_start < n_end; batch_start += batch_size) {
        /* Calculate actual batch size (may be smaller for last batch) */
        uint64_t actual_batch_size = batch_size;
        if (batch_start + batch_size > n_end) {
            actual_batch_size = n_end - batch_start;
        }

        /* Reset and process batch */
        batch_sieve_reset(bs, batch_start);
        bs->batch_size = actual_batch_size;
        batch_process(bs);

        /* Accumulate statistics */
        total_solved += bs->total_solved;
        total_mr_saved += bs->mr_tests_saved;
        total_mr_done += bs->mr_tests_done;
        batches_processed++;

        /* Check for counterexamples */
        if (bs->total_solved < actual_batch_size) {
            /* Verify unsolved n values */
            for (uint64_t idx = 0; idx < actual_batch_size; idx++) {
                if (!bs->solved[idx]) {
                    uint64_t n = batch_start + idx;
                    printf("\n*** POTENTIAL COUNTEREXAMPLE FOUND! ***\n");
                    printf("n = %s\n", fmt_num(n));
                    printf("N = 8n + 3 = %s\n", fmt_num(8*n + 3));
                    printf("Verifying with standard search...\n");

                    /* Double-check with standard algorithm */
                    uint64_t N = 8 * n + 3;
                    uint64_t a_max = isqrt64(N);
                    if ((a_max & 1) == 0) a_max--;

                    bool found = false;
                    for (uint64_t a = a_max; a >= 1; a -= 2) {
                        uint64_t a_sq = a * a;
                        if (a_sq > N - 4) continue;
                        uint64_t p = (N - a_sq) / 2;
                        if (p >= 2 && is_prime_64(p)) {
                            printf("VERIFIED: Solution exists! a=%llu, p=%llu\n",
                                   (unsigned long long)a, (unsigned long long)p);
                            found = true;
                            total_solved++;
                            break;
                        }
                    }

                    if (!found) {
                        printf("CONFIRMED: No solution exists!\n");
                        total_counterexamples++;
                    }
                }
            }
        }

        /* Progress reporting */
        double now;
#ifdef _OPENMP
        now = omp_get_wtime();
#else
        now = (double)clock() / CLOCKS_PER_SEC;
#endif
        double elapsed = now - global_start;

        if (elapsed - last_report_time >= PROGRESS_SECONDS) {
            double rate = total_solved / elapsed;
            double pct = 100.0 * total_solved / total;
            uint64_t remaining = total - total_solved;
            double eta_seconds = (rate > 0) ? remaining / rate : 0;

            printf("Batch %s/%s (%.1f%%), rate = %s n/sec, ETA: %s\n",
                   fmt_num(batches_processed), fmt_num(num_batches), pct,
                   fmt_num((uint64_t)rate), fmt_time(eta_seconds));
            fflush(stdout);

            last_report_time = elapsed;
        }
    }

    double global_end;
#ifdef _OPENMP
    global_end = omp_get_wtime();
#else
    global_end = (double)clock() / CLOCKS_PER_SEC;
#endif
    double global_elapsed = global_end - global_start;

    /* Print results */
    printf("\n");
    printf("==================================================================\n");
    printf("RESULTS\n");
    printf("==================================================================\n\n");

    printf("Total time:           %.2f seconds\n", global_elapsed);
    printf("Total throughput:     %s n/sec\n",
           fmt_num((uint64_t)(total / global_elapsed)));
    printf("Batches processed:    %s\n", fmt_num(batches_processed));
    printf("Counterexamples:      %s\n", fmt_num(total_counterexamples));

    printf("\nSieve Statistics:\n");
    uint64_t total_tests = total_mr_saved + total_mr_done;
    double save_rate = (total_tests > 0) ? 100.0 * total_mr_saved / total_tests : 0.0;
    printf("  MR tests saved:     %s (%.1f%%)\n", fmt_num(total_mr_saved), save_rate);
    printf("  MR tests performed: %s (%.1f%%)\n", fmt_num(total_mr_done), 100.0 - save_rate);

    /* Clean up */
    batch_sieve_destroy(bs);

    return (total_counterexamples > 0) ? 2 : 0;
}
