/*
 * GPU Counterexample Search: 8n + 3 = a^2 + 2p
 *
 * Metal GPU-accelerated version of the counterexample search.
 * Uses the GPU for parallel primality testing with CPU verification
 * of any potential counterexamples.
 *
 * Compile: make metal
 * Usage:   ./search_gpu [n_start] [n_end] [--batch-size N]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

/* Include shared headers */
#include "fmt.h"
#include "arith.h"
#include "prime.h"
#include "solve.h"

/* Metal host interface */
#include "metal_host.h"

/* Note: fj64_bases is already available via prime.h -> fj64_table.h */

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Progress reporting interval in seconds */
#define PROGRESS_SECONDS 5.0

/* Default search range */
#define DEFAULT_N_START 1000000000000ULL    /* 10^12 */
#define DEFAULT_N_END   1000010000000ULL    /* 10^12 + 10^7 */

/* Default batch size (can be overridden by --batch-size) */
#define DEFAULT_BATCH_SIZE 65536            /* 64K n values per GPU dispatch */

/* Threads per threadgroup (256 is typical for Metal) */
#define THREADS_PER_GROUP 256

/* ========================================================================== */
/* Time Utilities                                                             */
/* ========================================================================== */

/* Get wall clock time in seconds */
static double get_wall_time(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}

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
/* CPU Verification (Critical for correctness)                                */
/* ========================================================================== */

/**
 * Verify a potential counterexample on the CPU.
 *
 * GPU results MUST be verified on CPU before reporting as counterexamples.
 * This catches any GPU arithmetic errors or shader bugs.
 *
 * Returns: true if verified as counterexample, false if CPU found a solution
 */
static bool verify_counterexample_cpu(uint64_t n, uint64_t* out_a, uint64_t* out_p) {
    uint64_t p;
    uint64_t a = find_solution(n, &p);

    if (a > 0) {
        /* CPU found a solution - GPU was wrong */
        if (out_a) *out_a = a;
        if (out_p) *out_p = p;
        return false;
    }

    /* Both GPU and CPU agree: no solution exists */
    return true;
}

/* ========================================================================== */
/* GPU Search Driver                                                          */
/* ========================================================================== */

/**
 * Run GPU search over a range of n values with CPU verification.
 */
static uint64_t run_gpu_search(
    uint64_t n_start,
    uint64_t n_end,
    uint32_t batch_size,
    uint64_t* verified_counterexamples
) {
    uint64_t total = n_end - n_start;
    uint64_t processed = 0;
    uint64_t cpu_verified_counterexamples = 0;
    uint64_t gpu_false_positives = 0;

    /* Allocate buffers for batch processing */
    uint64_t* n_batch = (uint64_t*)malloc(batch_size * sizeof(uint64_t));
    GPUSearchResult* results = (GPUSearchResult*)malloc(batch_size * sizeof(GPUSearchResult));

    if (!n_batch || !results) {
        fprintf(stderr, "Error: Failed to allocate batch buffers\n");
        free(n_batch);
        free(results);
        return 0;
    }

    /* Get start time (wall clock) */
    double start_time = get_wall_time();
    double last_report_time = 0.0;

    printf("Processing %s n values in batches of %u...\n\n",
           fmt_num(total), batch_size);

    /* Process in batches */
    uint64_t current_n = n_start;
    while (current_n < n_end) {
        /* Fill batch */
        uint32_t batch_count = 0;
        while (batch_count < batch_size && current_n < n_end) {
            n_batch[batch_count++] = current_n++;
        }

        /* Execute GPU search */
        uint64_t batch_counterexamples = metal_search_batch(
            n_batch,
            batch_count,
            results,
            THREADS_PER_GROUP
        );

        processed += batch_count;
        (void)batch_counterexamples;  /* Used only for iteration below */

        /* Verify any potential counterexamples on CPU */
        for (uint32_t i = 0; i < batch_count; i++) {
            if (results[i].found == 0) {
                /* GPU reports no solution - verify on CPU */
                uint64_t cpu_a, cpu_p;
                bool is_genuine = verify_counterexample_cpu(results[i].n, &cpu_a, &cpu_p);

                if (is_genuine) {
                    /* Genuine counterexample confirmed! */
                    cpu_verified_counterexamples++;
                    printf("\n*** COUNTEREXAMPLE FOUND AND VERIFIED! ***\n");
                    printf("n = %s\n", fmt_num(results[i].n));
                    printf("N = 8n + 3 = %s\n", fmt_num(8 * results[i].n + 3));
                    printf("No valid (a, p) pair exists!\n");
                    printf("Verified by both GPU and CPU.\n\n");
                    fflush(stdout);
                } else {
                    /* GPU bug - CPU found a solution */
                    gpu_false_positives++;
                    printf("\nWARNING: GPU false positive for n = %s\n", fmt_num(results[i].n));
                    printf("  CPU found solution: a = %llu, p = %llu\n",
                           (unsigned long long)cpu_a, (unsigned long long)cpu_p);
                    printf("  Verification: 8n+3 = %llu, a^2+2p = %llu\n",
                           (unsigned long long)(8 * results[i].n + 3),
                           (unsigned long long)(cpu_a * cpu_a + 2 * cpu_p));
                    fflush(stdout);
                }
            }
        }

        /* Progress reporting */
        double now = get_wall_time();
        double elapsed = now - start_time;

        if (elapsed - last_report_time >= PROGRESS_SECONDS) {
            double rate = processed / elapsed;
            double pct = 100.0 * processed / total;
            uint64_t remaining = total - processed;
            double eta_seconds = (rate > 0) ? remaining / rate : 0;

            GPUStats stats;
            metal_get_stats(&stats);
            double gpu_rate = stats.total_n_processed / (stats.total_gpu_time_ms / 1000.0);

            printf("[GPU] n ~ %s (%.1f%%), rate = %s n/sec, GPU rate = %s n/sec, ETA: %s\n",
                   fmt_num(current_n), pct,
                   fmt_num((uint64_t)rate),
                   fmt_num((uint64_t)gpu_rate),
                   fmt_time(eta_seconds));
            fflush(stdout);

            last_report_time = elapsed;
        }
    }

    /* Cleanup */
    free(n_batch);
    free(results);

    /* Report false positive rate if any occurred */
    if (gpu_false_positives > 0) {
        printf("\nWARNING: %llu GPU false positives detected and corrected by CPU verification.\n",
               (unsigned long long)gpu_false_positives);
    }

    *verified_counterexamples = cpu_verified_counterexamples;
    return processed;
}

/* ========================================================================== */
/* Verification                                                               */
/* ========================================================================== */

/**
 * Verify the CPU algorithm against known solutions
 */
static bool verify_cpu_known_solutions(void) {
    printf("Verifying CPU algorithm...\n");

    struct { uint64_t n; uint64_t a; uint64_t p; } known[] = {
        {1, 1, 5},    /* 11 = 1 + 10 */
        {2, 3, 5},    /* 19 = 9 + 10 */
        {3, 1, 13},   /* 27 = 1 + 26 */
        {4, 5, 5}     /* 35 = 25 + 10 */
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

        uint64_t p_found;
        uint64_t found_a = find_solution(n, &p_found);

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

/**
 * Test GPU against CPU for a small range
 */
static bool verify_gpu_against_cpu(uint32_t test_count) {
    printf("Verifying GPU against CPU for %u values...\n", test_count);

    uint64_t* n_values = (uint64_t*)malloc(test_count * sizeof(uint64_t));
    GPUSearchResult* results = (GPUSearchResult*)malloc(test_count * sizeof(GPUSearchResult));

    if (!n_values || !results) {
        fprintf(stderr, "Error: Failed to allocate test buffers\n");
        free(n_values);
        free(results);
        return false;
    }

    /* Test range */
    for (uint32_t i = 0; i < test_count; i++) {
        n_values[i] = i + 1;
    }

    /* Run GPU search */
    metal_search_batch(n_values, test_count, results, THREADS_PER_GROUP);

    /* Compare with CPU */
    bool all_match = true;
    for (uint32_t i = 0; i < test_count; i++) {
        uint64_t n = n_values[i];
        uint64_t cpu_p;
        uint64_t cpu_a = find_solution(n, &cpu_p);

        bool cpu_found = (cpu_a > 0);
        bool gpu_found = (results[i].found != 0);

        if (cpu_found != gpu_found) {
            printf("  MISMATCH at n=%llu: CPU found=%d, GPU found=%d\n",
                   (unsigned long long)n, cpu_found, gpu_found);
            all_match = false;
        }
    }

    free(n_values);
    free(results);

    if (all_match) {
        printf("  All %u values match!\n", test_count);
    }

    return all_match;
}

/* ========================================================================== */
/* Argument Parsing                                                           */
/* ========================================================================== */

static uint64_t parse_number(const char* str) {
    char* endptr;
    double val = strtod(str, &endptr);
    if (*endptr != '\0') {
        return strtoull(str, NULL, 10);
    }
    return (uint64_t)val;
}

static void print_usage(const char* program) {
    printf("Usage: %s [n_start] [n_end] [options]\n", program);
    printf("\n");
    printf("GPU-accelerated search for counterexamples to: 8n + 3 = a^2 + 2p\n");
    printf("\n");
    printf("Uses Metal GPU compute for massive parallelization.\n");
    printf("All potential counterexamples are verified on CPU before reporting.\n");
    printf("\n");
    printf("Arguments:\n");
    printf("  n_start        Starting value of n (inclusive), default: 1e12\n");
    printf("  n_end          Ending value of n (exclusive), default: 1e12 + 1e7\n");
    printf("\n");
    printf("Options:\n");
    printf("  --batch-size N Number of n values per GPU dispatch (default: 65536)\n");
    printf("  --verify-only  Run verification tests only, don't search\n");
    printf("\n");
    printf("Numbers can be in scientific notation (e.g., 1e9, 2.5e6)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s                        Search [10^12, 10^12 + 10^7)\n", program);
    printf("  %s 1 1e6                  Search [1, 10^6)\n", program);
    printf("  %s 1e9 2e9 --batch-size 100000\n", program);
    printf("\n");
    printf("Exit codes:\n");
    printf("  0  Search completed, no counterexamples found\n");
    printf("  1  Error (no GPU, invalid arguments, verification failure)\n");
    printf("  2  Counterexample found and verified\n");
}

/* ========================================================================== */
/* Main                                                                       */
/* ========================================================================== */

int main(int argc, char** argv) {
    /* Disable stdout buffering for real-time output */
    setbuf(stdout, NULL);

    printf("==================================================================\n");
    printf("     GPU Counterexample Search: 8n + 3 = a^2 + 2p                 \n");
    printf("     Metal GPU Compute with CPU Verification                      \n");
    printf("==================================================================\n\n");

    uint64_t n_start = DEFAULT_N_START;
    uint64_t n_end = DEFAULT_N_END;
    uint32_t batch_size = DEFAULT_BATCH_SIZE;
    bool verify_only = false;

    /* Handle help flag */
    if (argc >= 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
        print_usage(argv[0]);
        return 0;
    }

    /* Parse arguments */
    int pos_count = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--batch-size") == 0 && i + 1 < argc) {
            batch_size = (uint32_t)parse_number(argv[++i]);
        } else if (strcmp(argv[i], "--verify-only") == 0) {
            verify_only = true;
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        } else {
            /* Positional argument */
            if (pos_count == 0) {
                n_start = parse_number(argv[i]);
            } else if (pos_count == 1) {
                n_end = parse_number(argv[i]);
            }
            pos_count++;
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

    /* Check for Metal availability */
    if (!metal_is_available()) {
        fprintf(stderr, "Error: Metal GPU not available\n");
        fprintf(stderr, "This program requires a Mac with Apple Silicon or AMD GPU.\n");
        return 1;
    }

    /* Initialize Metal */
    printf("Initializing Metal GPU...\n");
    if (!metal_init(fj64_bases)) {
        fprintf(stderr, "Error: Failed to initialize Metal\n");
        return 1;
    }

    printf("  Device: %s\n", metal_get_device_name());
    printf("  Compute units: %u\n", metal_get_compute_units());
    printf("  Recommended batch size: %u\n", metal_get_recommended_batch_size());
    printf("\n");

    /* Verify CPU algorithm */
    if (!verify_cpu_known_solutions()) {
        fprintf(stderr, "\nERROR: CPU verification failed!\n");
        metal_cleanup();
        return 1;
    }
    printf("\n");

    /* Verify GPU against CPU */
    if (!verify_gpu_against_cpu(1000)) {
        fprintf(stderr, "\nERROR: GPU verification failed!\n");
        metal_cleanup();
        return 1;
    }
    printf("\n");

    if (verify_only) {
        printf("Verification complete. Exiting.\n");
        metal_cleanup();
        return 0;
    }

    /* Print configuration */
    uint64_t total = n_end - n_start;
    printf("Configuration:\n");
    printf("  Range: n in [%s, %s)\n", fmt_num(n_start), fmt_num(n_end));
    printf("  Count: %s values\n", fmt_num(total));
    printf("  Batch size: %u\n", batch_size);
    printf("  Threads per group: %u\n", THREADS_PER_GROUP);
    printf("\n");

    /* Run GPU search */
    printf("Starting GPU search...\n\n");

    double global_start = get_wall_time();

    uint64_t verified_counterexamples = 0;
    uint64_t processed = run_gpu_search(n_start, n_end, batch_size, &verified_counterexamples);

    double global_end = get_wall_time();
    double global_elapsed = global_end - global_start;

    /* Get final statistics */
    GPUStats stats;
    metal_get_stats(&stats);

    /* Print results */
    printf("\n");
    printf("==================================================================\n");
    printf("RESULTS\n");
    printf("==================================================================\n\n");

    printf("Total time:           %.2f seconds\n", global_elapsed);
    printf("Values processed:     %s\n", fmt_num(processed));
    printf("Overall throughput:   %s n/sec\n", fmt_num((uint64_t)(total / global_elapsed)));
    printf("GPU time:             %.2f ms\n", stats.total_gpu_time_ms);
    printf("GPU throughput:       %s n/sec\n",
           fmt_num((uint64_t)(stats.total_n_processed / (stats.total_gpu_time_ms / 1000.0))));
    printf("Batches executed:     %llu\n", (unsigned long long)stats.total_batches);
    printf("GPU potential CEs:    %llu\n", (unsigned long long)stats.total_counterexamples);
    printf("Verified CEs:         %llu\n", (unsigned long long)verified_counterexamples);

    /* Cleanup */
    metal_cleanup();

    return (verified_counterexamples > 0) ? 2 : 0;
}
