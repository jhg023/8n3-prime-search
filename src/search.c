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
 * Parallelized with OpenMP for multi-core systems.
 *
 * Reference: Forisek & Jancina (2015), "Fast Primality Testing for
 * Integers That Fit into a Machine Word"
 *
 * Compile: make release
 * Usage:   ./search [n_start] [n_end] [--threads N]
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
#include "prime_sieve_fast.h"  /* Optimized: wheel30 + OpenMP (~20x faster) */

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Progress reporting interval in seconds */
#define PROGRESS_SECONDS 5.0

/* Default search range */
#define DEFAULT_N_START 1000000000000ULL    /* 10^12 */
#define DEFAULT_N_END   1000010000000ULL    /* 10^12 + 10^7 */

/* Maximum number of threads supported */
#define MAX_THREADS 256

/* ========================================================================== */
/* Time Formatting                                                            */
/* ========================================================================== */

/**
 * Format seconds into human-readable time string (e.g., "2h 15m 30s")
 * Uses a static buffer - not thread-safe, but only called from thread 0
 */
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
/* Per-Thread Statistics                                                      */
/* ========================================================================== */

typedef struct {
    uint64_t n_processed;      /* Number of n values processed */
    uint64_t total_checks;     /* Total a values checked */
    uint64_t counterexamples;  /* Counterexamples found by this thread */
    uint64_t sieve_hits;       /* Primes found via sieve lookup */
    uint64_t sieve_misses;     /* Primes found via Miller-Rabin */
    /* Padding to avoid false sharing (cache line is typically 64 bytes) */
    char padding[64 - 5 * sizeof(uint64_t)];
} ThreadStats;

static ThreadStats thread_stats[MAX_THREADS];

/* ========================================================================== */
/* Solution Finding (inlined for performance)                                 */
/* ========================================================================== */

/**
 * Check if candidate is filtered by trial division
 * Returns: 0 = composite (filtered), 1 = is small prime, 2 = needs Miller-Rabin
 */
static inline int trial_division_check_local(uint64_t candidate) {
    static const uint32_t PRIMES[] = {
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
        37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
        79, 83, 89, 97, 101, 103, 107, 109, 113, 127
    };

    /* Inline first 7 primes - catch ~65% of composites */
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    if (candidate % 11 == 0) return (candidate == 11) ? 1 : 0;
    if (candidate % 13 == 0) return (candidate == 13) ? 1 : 0;
    if (candidate % 17 == 0) return (candidate == 17) ? 1 : 0;
    if (candidate % 19 == 0) return (candidate == 19) ? 1 : 0;

    /* Check remaining primes 23-127 with 4x unrolled loop */
    int i = 7;
    while (i + 3 < 30) {
        if (candidate % PRIMES[i] == 0)
            return (candidate == PRIMES[i]) ? 1 : 0;
        if (candidate % PRIMES[i+1] == 0)
            return (candidate == PRIMES[i+1]) ? 1 : 0;
        if (candidate % PRIMES[i+2] == 0)
            return (candidate == PRIMES[i+2]) ? 1 : 0;
        if (candidate % PRIMES[i+3] == 0)
            return (candidate == PRIMES[i+3]) ? 1 : 0;
        i += 4;
    }
    while (i < 30) {
        if (candidate % PRIMES[i] == 0)
            return (candidate == PRIMES[i]) ? 1 : 0;
        i++;
    }
    return 2;
}

/**
 * Test if a candidate prime is actually prime
 */
static inline bool is_candidate_prime_local(uint64_t candidate) {
    int td = trial_division_check_local(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

/**
 * Test if a candidate prime is actually prime, using sieve when available.
 * Falls back to Miller-Rabin for candidates outside sieve range.
 */
static inline bool is_candidate_prime_with_sieve_local(uint64_t candidate,
                                                        const PrimeSieve *sieve,
                                                        int thread_id) {
    int td = trial_division_check_local(candidate);
    if (td == 0) return false;  /* Composite */
    if (td == 1) return true;   /* Small prime (3-127) */

    /* Candidate is > 127 and passed trial division */
    /* Try sieve lookup first if candidate is in range */
    if (sieve && sieve_in_range(sieve, candidate)) {
        thread_stats[thread_id].sieve_hits++;
        return sieve_is_prime(sieve, candidate);
    }

    /* Fall back to Miller-Rabin */
    thread_stats[thread_id].sieve_misses++;
    return is_prime_fj64_fast(candidate);
}

/**
 * Find a solution to 8n + 3 = a^2 + 2p
 * Returns the largest valid a, or 0 if no solution exists (counterexample).
 * Also updates thread-local statistics.
 */
static inline uint64_t find_solution_parallel(uint64_t n, int thread_id,
                                               const PrimeSieve *sieve) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Ensure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            thread_stats[thread_id].total_checks++;

            bool is_prime;
            if (sieve) {
                is_prime = is_candidate_prime_with_sieve_local(candidate, sieve, thread_id);
            } else {
                is_prime = is_candidate_prime_local(candidate);
            }

            if (is_prime) {
                thread_stats[thread_id].n_processed++;
                return a;
            }
        }

        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    thread_stats[thread_id].n_processed++;
    return 0;  /* Counterexample! */
}

/* ========================================================================== */
/* Parallel Search                                                            */
/* ========================================================================== */

/**
 * Run parallel search over a range of n values
 */
void run_search_parallel(uint64_t n_start, uint64_t n_end, int num_threads,
                         const PrimeSieve *sieve,
                         uint64_t *out_counterexamples) {
    uint64_t total_counterexamples = 0;
    uint64_t total = n_end - n_start;

    /* Initialize per-thread statistics */
    memset(thread_stats, 0, sizeof(thread_stats));

    /* Get start time */
    double start_time;
#ifdef _OPENMP
    start_time = omp_get_wtime();
#else
    start_time = (double)clock() / CLOCKS_PER_SEC;
#endif

    /* Progress reporting timing - shared across threads */
    volatile double last_report_time = 0.0;

    /* Early termination flag for counterexamples */
    volatile int found_counterexample = 0;

#ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel reduction(+:total_counterexamples)
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();
#else
        int tid = 0;
        int nthreads = 1;
#endif

        /* Calculate this thread's range */
        uint64_t chunk_size = (total + nthreads - 1) / nthreads;
        uint64_t my_start = n_start + (uint64_t)tid * chunk_size;
        uint64_t my_end = my_start + chunk_size;
        if (my_end > n_end) my_end = n_end;
        if (my_start >= n_end) my_start = my_end;  /* Empty range */

        uint64_t local_counterexamples = 0;
        uint64_t local_progress = 0;

        /* Process this thread's range */
        for (uint64_t n = my_start; n < my_end; n++) {
            /* Check for early termination */
            if (found_counterexample) break;

            uint64_t a = find_solution_parallel(n, tid, sieve);

            if (a == 0) {
                /* Counterexample found! */
                local_counterexamples++;
                found_counterexample = 1;  /* Signal all threads to stop */
#ifdef _OPENMP
                #pragma omp critical
#endif
                {
                    printf("\n*** COUNTEREXAMPLE FOUND! ***\n");
                    printf("n = %s (thread %d)\n", fmt_num(n), tid);
                    printf("N = 8n + 3 = %s\n", fmt_num(8*n + 3));
                    printf("No valid (a, p) pair exists!\n\n");
                    fflush(stdout);
                }
                break;  /* This thread stops immediately */
            }

            local_progress++;

            /* Progress reporting (any thread can report, with locking) */
            if ((local_progress & 0x3FFFF) == 0) {
                double now;
#ifdef _OPENMP
                now = omp_get_wtime();
#else
                now = (double)clock() / CLOCKS_PER_SEC;
#endif
                double elapsed = now - start_time;

                /* Only one thread reports at a time */
                if (elapsed - last_report_time >= PROGRESS_SECONDS) {
#ifdef _OPENMP
                    #pragma omp critical
#endif
                    {
                        /* Double-check timing inside critical section */
                        if (elapsed - last_report_time >= PROGRESS_SECONDS) {
                            /* Sum up all thread statistics */
                            uint64_t sum_processed = 0;
                            for (int t = 0; t < nthreads; t++) {
                                sum_processed += thread_stats[t].n_processed;
                            }

                            double rate = sum_processed / elapsed;
                            double pct = 100.0 * sum_processed / total;

                            /* Calculate ETA */
                            uint64_t remaining = total - sum_processed;
                            double eta_seconds = (rate > 0) ? remaining / rate : 0;

                            printf("[%d threads] n ~ %s (%.1f%%), rate = %s n/sec, ETA: %s\n",
                                   nthreads, fmt_num(n_start + sum_processed), pct,
                                   fmt_num((uint64_t)rate), fmt_time(eta_seconds));
                            fflush(stdout);

                            last_report_time = elapsed;
                        }
                    }
                }
            }
        }

        total_counterexamples += local_counterexamples;
    }

    *out_counterexamples = total_counterexamples;
}

/* ========================================================================== */
/* Verification                                                               */
/* ========================================================================== */

/**
 * Verify the algorithm against known solutions
 */
bool verify_known_solutions(const PrimeSieve *sieve) {
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

        uint64_t found_a = find_solution_parallel(n, 0, sieve);

        printf("  n=%llu: N=%llu, given (%llu,%llu), found a=%llu ... ",
               (unsigned long long)n, (unsigned long long)N,
               (unsigned long long)expected_a, (unsigned long long)expected_p,
               (unsigned long long)found_a);

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
    printf("Usage: %s [n_start] [n_end] [options]\n", program);
    printf("\n");
    printf("Search for counterexamples to: 8n + 3 = a^2 + 2p\n");
    printf("\n");
    printf("For each n in [n_start, n_end), attempts to find odd a and prime p\n");
    printf("such that 8n + 3 = a^2 + 2p. Reports any n for which no solution exists.\n");
    printf("\n");
    printf("Arguments:\n");
    printf("  n_start              Starting value of n (inclusive), default: 1e12\n");
    printf("  n_end                Ending value of n (exclusive), default: 1e12 + 1e7\n");
    printf("  --threads N          Number of threads to use (default: all cores)\n");
    printf("  --sieve-threshold T  Pre-compute prime sieve up to T for O(1) lookups\n");
    printf("                       Recommended values: 1e7 (1MB), 1e8 (12MB), 1e9 (125MB)\n");
    printf("\n");
    printf("Numbers can be in scientific notation (e.g., 1e9, 2.5e6)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s                        Search [10^12, 10^12 + 10^7) with all cores\n", program);
    printf("  %s 1 1e6                  Search [1, 10^6)\n", program);
    printf("  %s 1e9 2e9 --threads 4    Search [10^9, 2*10^9) with 4 threads\n", program);
    printf("  %s 1e12 1.001e12 --sieve-threshold 1e8  Use 12MB prime sieve\n", program);
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
    printf("     Optimized with FJ64_262K primality test (OpenMP parallel)    \n");
    printf("==================================================================\n\n");

    uint64_t n_start = DEFAULT_N_START;
    uint64_t n_end = DEFAULT_N_END;
    int num_threads = 0;  /* 0 = auto-detect */
    uint64_t sieve_threshold = 0;  /* 0 = no sieve */

    /* Handle help flag */
    if (argc >= 2 && (strcmp(argv[1], "-h") == 0 ||
                      strcmp(argv[1], "--help") == 0)) {
        print_usage(argv[0]);
        return 0;
    }

    /* Parse arguments */
    int arg_idx = 1;
    while (arg_idx < argc) {
        if (strcmp(argv[arg_idx], "--threads") == 0 && arg_idx + 1 < argc) {
            num_threads = atoi(argv[arg_idx + 1]);
            arg_idx += 2;
        } else if (strcmp(argv[arg_idx], "--sieve-threshold") == 0 && arg_idx + 1 < argc) {
            sieve_threshold = parse_number(argv[arg_idx + 1]);
            arg_idx += 2;
        } else if (argv[arg_idx][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[arg_idx]);
            print_usage(argv[0]);
            return 1;
        } else {
            /* Positional argument */
            if (n_start == DEFAULT_N_START && arg_idx == 1) {
                n_start = parse_number(argv[arg_idx]);
                n_end = n_start + 10000000;  /* Default range if only start given */
            } else if (arg_idx == 2 || (arg_idx == 1 && n_start != DEFAULT_N_START)) {
                n_end = parse_number(argv[arg_idx]);
            }
            arg_idx++;
        }
    }

    /* Re-parse positional args more carefully */
    arg_idx = 1;
    int pos_count = 0;
    while (arg_idx < argc) {
        if (strcmp(argv[arg_idx], "--threads") == 0 ||
            strcmp(argv[arg_idx], "--sieve-threshold") == 0) {
            arg_idx += 2;
            continue;
        }
        if (argv[arg_idx][0] == '-') {
            arg_idx++;
            continue;
        }
        if (pos_count == 0) {
            n_start = parse_number(argv[arg_idx]);
        } else if (pos_count == 1) {
            n_end = parse_number(argv[arg_idx]);
        }
        pos_count++;
        arg_idx++;
    }

    /* Set default end if only start was given */
    if (pos_count == 1) {
        n_end = n_start + 10000000;
    }

    if (n_start >= n_end) {
        fprintf(stderr, "Error: n_start must be less than n_end\n");
        return 1;
    }

    /* Determine number of threads */
#ifdef _OPENMP
    if (num_threads <= 0) {
        num_threads = omp_get_max_threads();
    }
    if (num_threads > MAX_THREADS) {
        num_threads = MAX_THREADS;
    }
#else
    num_threads = 1;
    if (num_threads > 1) {
        fprintf(stderr, "Warning: OpenMP not available, running single-threaded\n");
    }
#endif

    uint64_t total = n_end - n_start;

    /* Create prime sieve if requested */
    PrimeSieve *sieve = NULL;
    if (sieve_threshold > 0) {
        printf("Creating prime sieve up to %s...\n", fmt_num(sieve_threshold));
        double sieve_start;
#ifdef _OPENMP
        sieve_start = omp_get_wtime();
#else
        sieve_start = (double)clock() / CLOCKS_PER_SEC;
#endif
        sieve = sieve_create(sieve_threshold);
        double sieve_end;
#ifdef _OPENMP
        sieve_end = omp_get_wtime();
#else
        sieve_end = (double)clock() / CLOCKS_PER_SEC;
#endif
        if (!sieve) {
            fprintf(stderr, "Error: Failed to allocate prime sieve\n");
            return 1;
        }
        printf("  Sieve created in %.2f seconds\n", sieve_end - sieve_start);
        printf("  Memory usage: %s\n", sieve_memory_str(sieve));
        printf("  Primes found: %s\n", fmt_num(sieve_prime_count(sieve)));
        printf("\n");
    }

    printf("Configuration:\n");
    printf("  Range: n in [%s, %s)\n", fmt_num(n_start), fmt_num(n_end));
    printf("  Count: %s values\n", fmt_num(total));
    printf("  Threads: %d\n", num_threads);
    if (sieve) {
        printf("  Primality test: Sieve lookup (up to %s) + FJ64_262K\n",
               fmt_num(sieve_threshold));
    } else {
        printf("  Primality test: FJ64_262K (2 Miller-Rabin tests)\n");
    }
    printf("\n");

    /* Verify algorithm correctness */
    if (!verify_known_solutions(sieve)) {
        fprintf(stderr, "\nERROR: Verification failed!\n");
        if (sieve) sieve_destroy(sieve);
        return 1;
    }
    printf("\n");

    /* Run search */
    printf("Starting parallel search...\n\n");

    double global_start;
#ifdef _OPENMP
    global_start = omp_get_wtime();
#else
    global_start = (double)clock() / CLOCKS_PER_SEC;
#endif

    uint64_t total_counterexamples = 0;
    run_search_parallel(n_start, n_end, num_threads, sieve, &total_counterexamples);

    double global_end;
#ifdef _OPENMP
    global_end = omp_get_wtime();
#else
    global_end = (double)clock() / CLOCKS_PER_SEC;
#endif
    double global_elapsed = global_end - global_start;

    /* Sum final statistics */
    uint64_t stat_n = 0, stat_checks = 0;
    uint64_t stat_sieve_hits = 0, stat_sieve_misses = 0;
    for (int t = 0; t < num_threads; t++) {
        stat_n += thread_stats[t].n_processed;
        stat_checks += thread_stats[t].total_checks;
        stat_sieve_hits += thread_stats[t].sieve_hits;
        stat_sieve_misses += thread_stats[t].sieve_misses;
    }
    double avg_checks = (stat_n > 0) ? (double)stat_checks / stat_n : 0.0;

    /* Print results */
    printf("\n");
    printf("==================================================================\n");
    printf("RESULTS\n");
    printf("==================================================================\n\n");

    printf("Total time:           %.2f seconds\n", global_elapsed);
    printf("Threads used:         %d\n", num_threads);
    printf("Total throughput:     %s n/sec\n",
           fmt_num((uint64_t)(total / global_elapsed)));
    printf("Per-thread avg:       %s n/sec\n",
           fmt_num((uint64_t)(total / global_elapsed / num_threads)));
    printf("Counterexamples:      %s\n", fmt_num(total_counterexamples));
    printf("Avg checks per n:     %.2f\n", avg_checks);
    printf("Total a's checked:    %s\n", fmt_num(stat_checks));

    /* Sieve statistics */
    if (sieve) {
        uint64_t total_primes_found = stat_sieve_hits + stat_sieve_misses;
        double hit_rate = (total_primes_found > 0)
            ? 100.0 * stat_sieve_hits / total_primes_found : 0.0;
        printf("\nSieve Statistics:\n");
        printf("  Sieve hits:         %s (%.1f%%)\n", fmt_num(stat_sieve_hits), hit_rate);
        printf("  Miller-Rabin tests: %s (%.1f%%)\n", fmt_num(stat_sieve_misses), 100.0 - hit_rate);
        printf("  Sieve threshold:    %s\n", fmt_num(sieve_threshold));
    }

    /* Clean up */
    if (sieve) {
        sieve_destroy(sieve);
    }

    return (total_counterexamples > 0) ? 2 : 0;
}
