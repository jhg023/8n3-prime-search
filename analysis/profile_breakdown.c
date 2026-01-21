/*
 * Profile where time is spent in the search
 * Measures: iteration overhead, trial division, Miller-Rabin
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "arith.h"
#include "prime.h"

/* Use clock_gettime for high-resolution timing on macOS */
static inline uint64_t get_ns(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
}

/* Counters */
static uint64_t cnt_candidates = 0;
static uint64_t cnt_trial_filtered = 0;
static uint64_t cnt_mr_tested = 0;
static uint64_t cnt_primes_found = 0;

static uint64_t ns_iteration = 0;
static uint64_t ns_trial = 0;
static uint64_t ns_mr = 0;

/* Modified primality test with profiling */
static inline bool is_candidate_prime_profiled(uint64_t candidate) {
    uint64_t t0 = get_ns();

    /* Trial division */
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            uint64_t t1 = get_ns();
            ns_trial += t1 - t0;
            cnt_trial_filtered++;
            return (candidate == TRIAL_PRIMES[i]);
        }
    }

    uint64_t t1 = get_ns();
    ns_trial += t1 - t0;

    /* Small primes that passed trial */
    if (candidate <= 127) {
        cnt_primes_found++;
        return true;
    }

    /* Miller-Rabin */
    cnt_mr_tested++;
    bool is_prime = is_prime_fj64_fast(candidate);
    uint64_t t2 = get_ns();
    ns_mr += t2 - t1;

    if (is_prime) cnt_primes_found++;
    return is_prime;
}

/* Find solution with profiling */
static uint64_t find_solution_profiled(uint64_t n) {
    uint64_t t_start = get_ns();

    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;
        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;
            cnt_candidates++;

            uint64_t t_iter = get_ns();
            ns_iteration += t_iter - t_start;

            if (is_candidate_prime_profiled(candidate)) {
                return a;
            }

            t_start = get_ns();
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

int main(void) {
    printf("Profiling Search Bottlenecks\n");
    printf("============================\n\n");

    uint64_t test_ns[] = {1000000000000ULL, 2000000000000000000ULL};
    const char* labels[] = {"10^12", "2e18"};
    uint64_t counts[] = {100000, 10000};

    for (int t = 0; t < 2; t++) {
        uint64_t n_start = test_ns[t];
        uint64_t count = counts[t];

        /* Reset counters */
        cnt_candidates = cnt_trial_filtered = cnt_mr_tested = cnt_primes_found = 0;
        ns_iteration = ns_trial = ns_mr = 0;

        printf("=== n ~ %s (testing %llu values) ===\n\n", labels[t], (unsigned long long)count);

        clock_t wall_start = clock();

        for (uint64_t i = 0; i < count; i++) {
            find_solution_profiled(n_start + i);
        }

        clock_t wall_end = clock();
        double wall_time = (double)(wall_end - wall_start) / CLOCKS_PER_SEC;

        uint64_t total_ns = ns_iteration + ns_trial + ns_mr;

        printf("Counts:\n");
        printf("  n values tested:      %llu\n", (unsigned long long)count);
        printf("  Candidates checked:   %llu (%.2f per n)\n",
               (unsigned long long)cnt_candidates, (double)cnt_candidates / count);
        printf("  Trial-div filtered:   %llu (%.1f%% of candidates)\n",
               (unsigned long long)cnt_trial_filtered,
               100.0 * cnt_trial_filtered / cnt_candidates);
        printf("  Miller-Rabin tested:  %llu (%.1f%% of candidates)\n",
               (unsigned long long)cnt_mr_tested,
               100.0 * cnt_mr_tested / cnt_candidates);
        printf("  Primes found:         %llu\n", (unsigned long long)cnt_primes_found);

        printf("\nCycle breakdown:\n");
        printf("  Iteration overhead:   %.1f%% (%llu ns)\n",
               100.0 * ns_iteration / total_ns, (unsigned long long)ns_iteration);
        printf("  Trial division:       %.1f%% (%llu ns)\n",
               100.0 * ns_trial / total_ns, (unsigned long long)ns_trial);
        printf("  Miller-Rabin:         %.1f%% (%llu ns)\n",
               100.0 * ns_mr / total_ns, (unsigned long long)ns_mr);

        printf("\nThroughput: %.0f n/sec (wall time: %.2fs)\n\n", count / wall_time, wall_time);
    }

    return 0;
}
