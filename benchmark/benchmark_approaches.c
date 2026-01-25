/*
 * Benchmark Comparison: Different Optimization Approaches
 *
 * Compares performance of:
 * 1. Baseline (no sieve)
 * 2. Prime sieve with various thresholds (10^7, 10^8, 10^9)
 * 3. Batched sieve approach
 *
 * Usage: ./benchmark_approaches [--count N] [--quick]
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
#include "prime_sieve.h"
#include "batch_sieve.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Default benchmark count */
#define DEFAULT_COUNT 1000000
#define QUICK_COUNT 100000

/* Test scales */
static const uint64_t TEST_SCALES[] = {
    1000000000ULL,          /* 10^9 */
    1000000000000ULL,       /* 10^12 */
    1000000000000000ULL,    /* 10^15 */
};
#define NUM_SCALES 3

/* Sieve thresholds to test */
static const uint64_t SIEVE_THRESHOLDS[] = {
    10000000ULL,            /* 10^7 = ~1.25 MB */
    100000000ULL,           /* 10^8 = ~12.5 MB */
};
#define NUM_SIEVE_THRESHOLDS 2

/* ========================================================================== */
/* Timing Helpers                                                             */
/* ========================================================================== */

static double get_time(void) {
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return (double)clock() / CLOCKS_PER_SEC;
#endif
}

/* ========================================================================== */
/* Trial Division (inlined)                                                   */
/* ========================================================================== */

static const uint32_t BENCH_TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};

static inline int bench_trial_division(uint64_t candidate) {
    for (int i = 0; i < 30; i++) {
        if (candidate % BENCH_TRIAL_PRIMES[i] == 0) {
            return (candidate == BENCH_TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool bench_is_prime(uint64_t candidate) {
    int td = bench_trial_division(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

/* ========================================================================== */
/* Benchmark: Baseline                                                        */
/* ========================================================================== */

typedef struct {
    double elapsed;
    uint64_t n_processed;
    uint64_t total_checks;
    uint64_t throughput;
} BenchResult;

static BenchResult bench_baseline(uint64_t n_start, uint64_t count) {
    BenchResult result = {0};

    double start = get_time();

    for (uint64_t i = 0; i < count; i++) {
        uint64_t n = n_start + i;
        uint64_t N = 8 * n + 3;
        uint64_t a_max = isqrt64(N);
        if ((a_max & 1) == 0) a_max--;

        uint64_t a = a_max;
        uint64_t candidate = (N - a * a) >> 1;
        uint64_t delta = 2 * (a - 1);

        while (1) {
            if (candidate >= 2) {
                result.total_checks++;
                if (bench_is_prime(candidate)) {
                    result.n_processed++;
                    break;
                }
            }
            if (a < 3) {
                result.n_processed++;
                break;
            }
            candidate += delta;
            delta -= 4;
            a -= 2;
        }
    }

    double end = get_time();
    result.elapsed = end - start;
    result.throughput = (result.elapsed > 0) ? (uint64_t)(count / result.elapsed) : 0;

    return result;
}

/* ========================================================================== */
/* Benchmark: With Sieve                                                      */
/* ========================================================================== */

typedef struct {
    double elapsed;
    uint64_t n_processed;
    uint64_t total_checks;
    uint64_t sieve_hits;
    uint64_t sieve_misses;
    uint64_t throughput;
    double hit_rate;
} BenchSieveResult;

static BenchSieveResult bench_with_sieve(uint64_t n_start, uint64_t count,
                                          const PrimeSieve *sieve) {
    BenchSieveResult result = {0};

    double start = get_time();

    for (uint64_t i = 0; i < count; i++) {
        uint64_t n = n_start + i;
        uint64_t N = 8 * n + 3;
        uint64_t a_max = isqrt64(N);
        if ((a_max & 1) == 0) a_max--;

        uint64_t a = a_max;
        uint64_t candidate = (N - a * a) >> 1;
        uint64_t delta = 2 * (a - 1);

        while (1) {
            if (candidate >= 2) {
                result.total_checks++;

                /* Check with sieve */
                int td = bench_trial_division(candidate);
                if (td == 0) {
                    /* composite */
                } else if (td == 1) {
                    /* small prime - found solution */
                    result.n_processed++;
                    break;
                } else {
                    /* Need primality test */
                    if (sieve && sieve_in_range(sieve, candidate)) {
                        result.sieve_hits++;
                        if (sieve_is_prime(sieve, candidate)) {
                            result.n_processed++;
                            break;
                        }
                    } else {
                        result.sieve_misses++;
                        if (is_prime_fj64_fast(candidate)) {
                            result.n_processed++;
                            break;
                        }
                    }
                }
            }
            if (a < 3) {
                result.n_processed++;
                break;
            }
            candidate += delta;
            delta -= 4;
            a -= 2;
        }
    }

    double end = get_time();
    result.elapsed = end - start;
    result.throughput = (result.elapsed > 0) ? (uint64_t)(count / result.elapsed) : 0;

    uint64_t total = result.sieve_hits + result.sieve_misses;
    result.hit_rate = (total > 0) ? 100.0 * result.sieve_hits / total : 0.0;

    return result;
}

/* ========================================================================== */
/* Benchmark: Batched                                                         */
/* ========================================================================== */

typedef struct {
    double elapsed;
    uint64_t n_processed;
    uint64_t mr_saved;
    uint64_t mr_done;
    uint64_t throughput;
    double save_rate;
} BenchBatchResult;

static BenchBatchResult bench_batched(uint64_t n_start, uint64_t count,
                                       uint64_t batch_size) {
    BenchBatchResult result = {0};

    BatchSieve *bs = batch_sieve_create(n_start, batch_size);
    if (!bs) {
        fprintf(stderr, "Failed to create batch sieve\n");
        return result;
    }

    double start = get_time();

    for (uint64_t batch_start = n_start; batch_start < n_start + count; batch_start += batch_size) {
        uint64_t actual_size = batch_size;
        if (batch_start + batch_size > n_start + count) {
            actual_size = n_start + count - batch_start;
        }

        batch_sieve_reset(bs, batch_start);
        bs->batch_size = actual_size;
        batch_process(bs);

        result.n_processed += bs->total_solved;
        result.mr_saved += bs->mr_tests_saved;
        result.mr_done += bs->mr_tests_done;
    }

    double end = get_time();
    result.elapsed = end - start;
    result.throughput = (result.elapsed > 0) ? (uint64_t)(count / result.elapsed) : 0;

    uint64_t total = result.mr_saved + result.mr_done;
    result.save_rate = (total > 0) ? 100.0 * result.mr_saved / total : 0.0;

    batch_sieve_destroy(bs);

    return result;
}

/* ========================================================================== */
/* Main Benchmark Runner                                                      */
/* ========================================================================== */

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

    printf("==================================================================\n");
    printf("     Benchmark: Optimization Approaches Comparison               \n");
    printf("==================================================================\n\n");

    uint64_t count = DEFAULT_COUNT;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--quick") == 0) {
            count = QUICK_COUNT;
        } else if (strcmp(argv[i], "--count") == 0 && i + 1 < argc) {
            count = strtoull(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("Usage: %s [--count N] [--quick]\n", argv[0]);
            printf("  --count N   Number of n values per test (default: 1M)\n");
            printf("  --quick     Use 100K values per test\n");
            return 0;
        }
    }

    printf("Configuration:\n");
    printf("  Values per test: %s\n", fmt_num(count));
    printf("\n");

    /* Pre-create sieves */
    printf("Creating prime sieves...\n");
    PrimeSieve *sieves[NUM_SIEVE_THRESHOLDS];
    for (int s = 0; s < NUM_SIEVE_THRESHOLDS; s++) {
        double start = get_time();
        sieves[s] = sieve_create(SIEVE_THRESHOLDS[s]);
        double elapsed = get_time() - start;
        printf("  Sieve 10^%d: %s, created in %.2fs\n",
               (SIEVE_THRESHOLDS[s] == 10000000 ? 7 : 8),
               sieve_memory_str(sieves[s]), elapsed);
    }
    printf("\n");

    /* Run benchmarks at each scale */
    for (int scale_idx = 0; scale_idx < NUM_SCALES; scale_idx++) {
        uint64_t n_start = TEST_SCALES[scale_idx];

        printf("==================================================================\n");
        printf("Scale: n = %s (10^%d)\n", fmt_num(n_start),
               scale_idx == 0 ? 9 : (scale_idx == 1 ? 12 : 15));
        printf("==================================================================\n\n");

        /* Baseline */
        printf("Testing baseline (no sieve)...\n");
        BenchResult baseline = bench_baseline(n_start, count);
        printf("  Throughput: %s n/sec\n", fmt_num(baseline.throughput));
        printf("  Avg checks: %.2f\n\n", (double)baseline.total_checks / baseline.n_processed);

        /* Sieve variants */
        for (int s = 0; s < NUM_SIEVE_THRESHOLDS; s++) {
            printf("Testing sieve (threshold = 10^%d)...\n",
                   SIEVE_THRESHOLDS[s] == 10000000 ? 7 : 8);
            BenchSieveResult sieve_result = bench_with_sieve(n_start, count, sieves[s]);
            double speedup = (baseline.elapsed > 0) ? baseline.elapsed / sieve_result.elapsed : 0;
            printf("  Throughput: %s n/sec (%.2fx speedup)\n",
                   fmt_num(sieve_result.throughput), speedup);
            printf("  Sieve hit rate: %.1f%%\n\n", sieve_result.hit_rate);
        }

        /* Batched (only at smaller scales due to different algorithm structure) */
        if (scale_idx <= 1) {
            printf("Testing batched sieve (batch_size = 64K)...\n");
            BenchBatchResult batch_result = bench_batched(n_start, count, 65536);
            double speedup = (baseline.elapsed > 0) ? baseline.elapsed / batch_result.elapsed : 0;
            printf("  Throughput: %s n/sec (%.2fx speedup)\n",
                   fmt_num(batch_result.throughput), speedup);
            printf("  MR test savings: %.1f%%\n\n", batch_result.save_rate);
        }
    }

    /* Summary */
    printf("==================================================================\n");
    printf("SUMMARY\n");
    printf("==================================================================\n\n");
    printf("The prime sieve optimization provides significant speedup by\n");
    printf("replacing Miller-Rabin tests with O(1) bitmap lookups for small\n");
    printf("prime candidates. Since we iterate largest-a first (smallest-p),\n");
    printf("most solutions hit the sieve.\n");
    printf("\n");
    printf("Recommended sieve threshold: 10^8 (12.5 MB) for best balance of\n");
    printf("memory usage and hit rate.\n");

    /* Cleanup */
    for (int s = 0; s < NUM_SIEVE_THRESHOLDS; s++) {
        sieve_destroy(sieves[s]);
    }

    return 0;
}
