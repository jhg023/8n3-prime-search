/*
 * Benchmark Suite for Counterexample Search: 8n + 3 = a^2 + 2p
 *
 * Tests throughput at various scales from 10^6 to ~2^61 for consistent
 * comparison across code changes.
 *
 * Compile: make benchmark
 * Usage:   ./benchmark_suite [--quick] [--count N]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#include "fmt.h"
#include "arith.h"
#include "solve.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

#define DEFAULT_COUNT   10000000  /* 10M iterations per scale */
#define QUICK_COUNT     1000000   /* 1M for quick mode */
#define WARMUP_COUNT    100000    /* 100k warmup iterations */

/* Test scales from 10^6 to ~2^61 */
static const struct {
    uint64_t n_start;
    const char* label;
    int bits;  /* Approximate bits of N = 8n+3 */
} SCALES[] = {
    {1000000ULL,                   "10^6",   23},
    {1000000000ULL,                "10^9",   33},
    {1000000000000ULL,             "10^12",  43},
    {1000000000000000ULL,          "10^15",  53},
    {100000000000000000ULL,        "10^17",  60},
    {1000000000000000000ULL,       "10^18",  63},
    {2000000000000000000ULL,       "2e18",   64},  /* ~2^61 for n, ~2^64 for N */
};
#define NUM_SCALES (sizeof(SCALES) / sizeof(SCALES[0]))

/* ========================================================================== */
/* Timing                                                                     */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Benchmark Core                                                             */
/* ========================================================================== */

typedef struct {
    uint64_t n_start;
    uint64_t count;
    double elapsed_sec;
    double n_per_sec;
    double avg_checks;  /* Average number of a values checked per n */
} BenchResult;

BenchResult run_benchmark(uint64_t n_start, uint64_t count) {
    BenchResult result = {0};
    result.n_start = n_start;
    result.count = count;

    /* Warmup */
    for (uint64_t i = 0; i < WARMUP_COUNT && i < count; i++) {
        uint64_t p;
        volatile uint64_t a = find_solution(n_start + i, &p);
        (void)a;
    }

    /* Timed run - also count how many a values we check */
    uint64_t total_checks = 0;
    double start = get_time();

    for (uint64_t i = 0; i < count; i++) {
        uint64_t n = n_start + i;
        uint64_t N = 8 * n + 3;
        uint64_t a_max = isqrt64(N);
        if ((a_max & 1) == 0) a_max--;

        uint64_t p;
        uint64_t a_found = find_solution(n, &p);

        /*
         * Count checks: we iterate from a_max down to a_found
         * Number of odd values from a_found to a_max = (a_max - a_found) / 2 + 1
         */
        if (a_found > 0) {
            total_checks += (a_max - a_found) / 2 + 1;
        }
    }

    double end = get_time();

    result.elapsed_sec = end - start;
    result.n_per_sec = count / result.elapsed_sec;
    result.avg_checks = (double)total_checks / count;

    return result;
}

/* ========================================================================== */
/* Main                                                                       */
/* ========================================================================== */

void print_usage(const char* program) {
    printf("Usage: %s [OPTIONS]\n\n", program);
    printf("Options:\n");
    printf("  --quick       Run with 1M iterations (faster)\n");
    printf("  --count N     Set iterations per scale (default: 10M)\n");
    printf("  -h, --help    Show this help message\n");
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

    uint64_t count = DEFAULT_COUNT;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--quick") == 0) {
            count = QUICK_COUNT;
        } else if (strcmp(argv[i], "--count") == 0 && i + 1 < argc) {
            count = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            print_usage(argv[0]);
            return 0;
        } else {
            fprintf(stderr, "Unknown option: %s\n", argv[i]);
            print_usage(argv[0]);
            return 1;
        }
    }

    /* Print header */
    printf("Benchmark: 8n + 3 = a^2 + 2p\n");
    printf("Iterations per scale: %s\n\n", fmt_num(count));

    printf("%-8s  %6s  %15s  %12s  %8s\n",
           "Scale", "Bits", "Rate (n/sec)", "Avg checks", "Time (s)");
    printf("--------------------------------------------------------------\n");

    /* Run benchmarks */
    for (size_t i = 0; i < NUM_SCALES; i++) {
        BenchResult res = run_benchmark(SCALES[i].n_start, count);

        printf("%-8s  %6d  %15s  %12.2f  %8.2f\n",
               SCALES[i].label,
               SCALES[i].bits,
               fmt_num((uint64_t)res.n_per_sec),
               res.avg_checks,
               res.elapsed_sec);
    }

    printf("--------------------------------------------------------------\n");

    return 0;
}
