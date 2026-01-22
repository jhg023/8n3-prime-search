/*
 * Per-line cycle profiling using rdtsc (or clock_gettime on ARM)
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>

#include "arith.h"
#include "fj64_table.h"

/* High-resolution cycle counter */
#if defined(__x86_64__)
static inline uint64_t rdtsc(void) {
    uint32_t lo, hi;
    __asm__ volatile ("rdtsc" : "=a"(lo), "=d"(hi));
    return ((uint64_t)hi << 32) | lo;
}
#else
/* ARM: use clock_gettime */
static inline uint64_t rdtsc(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
}
#endif

/* Profiling counters */
static uint64_t cycles_isqrt = 0;
static uint64_t cycles_loop_overhead = 0;
static uint64_t cycles_a_squared = 0;
static uint64_t cycles_candidate_calc = 0;
static uint64_t cycles_trial_div = 0;
static uint64_t cycles_montgomery_setup = 0;
static uint64_t cycles_montgomery_exp = 0;
static uint64_t cycles_montgomery_witness = 0;
static uint64_t cycles_total = 0;

static uint64_t count_candidates = 0;
static uint64_t count_trial_filtered = 0;
static uint64_t count_mr_tests = 0;

/* Trial primes */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_TRIAL_PRIMES 30

/* Montgomery constants */
#define MONTGOMERY_SAFE_THRESHOLD (1ULL << 63)

static inline uint64_t montgomery_inverse(uint64_t n) {
    uint64_t x = n;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    return -x;
}

static inline uint64_t montgomery_reduce(__uint128_t t, uint64_t n, uint64_t n_inv) {
    uint64_t m = (uint64_t)t * n_inv;
    __uint128_t mn = (__uint128_t)m * n;
    uint64_t u = (t + mn) >> 64;
    return (u >= n) ? (u - n) : u;
}

static inline uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t n, uint64_t n_inv) {
    return montgomery_reduce((__uint128_t)a * b, n, n_inv);
}

static inline uint64_t montgomery_r_squared(uint64_t n) {
    uint64_t r = (((__uint128_t)1 << 64) % n);
    return ((__uint128_t)r * r) % n;
}

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    return ((x >> 32) ^ x) & 262143;
}

/* Profiled Miller-Rabin */
static inline bool mr_witness_profiled(uint64_t n, uint64_t a) {
    uint64_t t0, t1;

    if (a >= n) a %= n;
    if (a == 0) return true;

    uint64_t d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;

    /* Montgomery setup */
    t0 = rdtsc();
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    uint64_t one_m = montgomery_reduce(r_sq, n, n_inv);
    uint64_t neg_one_m = n - one_m;
    uint64_t a_m = montgomery_reduce((__uint128_t)a * r_sq, n, n_inv);
    t1 = rdtsc();
    cycles_montgomery_setup += t1 - t0;

    /* Exponentiation */
    t0 = rdtsc();
    uint64_t x_m = one_m;
    uint64_t base_m = a_m;
    uint64_t exp = d;
    while (exp > 0) {
        if (exp & 1) {
            x_m = montgomery_mul(x_m, base_m, n, n_inv);
        }
        base_m = montgomery_mul(base_m, base_m, n, n_inv);
        exp >>= 1;
    }
    t1 = rdtsc();
    cycles_montgomery_exp += t1 - t0;

    /* Witness loop */
    t0 = rdtsc();
    if (x_m == one_m || x_m == neg_one_m) {
        t1 = rdtsc();
        cycles_montgomery_witness += t1 - t0;
        return true;
    }

    for (int i = 1; i < r; i++) {
        x_m = montgomery_mul(x_m, x_m, n, n_inv);
        if (x_m == neg_one_m) {
            t1 = rdtsc();
            cycles_montgomery_witness += t1 - t0;
            return true;
        }
        if (x_m == one_m) {
            t1 = rdtsc();
            cycles_montgomery_witness += t1 - t0;
            return false;
        }
    }
    t1 = rdtsc();
    cycles_montgomery_witness += t1 - t0;
    return false;
}

/* Profiled primality test */
static inline bool is_prime_profiled(uint64_t candidate) {
    uint64_t t0, t1;

    /* Trial division */
    t0 = rdtsc();
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            t1 = rdtsc();
            cycles_trial_div += t1 - t0;
            count_trial_filtered++;
            return (candidate == TRIAL_PRIMES[i]);
        }
    }
    t1 = rdtsc();
    cycles_trial_div += t1 - t0;

    if (candidate <= 127) return true;

    /* Miller-Rabin */
    count_mr_tests++;
    if (!mr_witness_profiled(candidate, 2))
        return false;
    return mr_witness_profiled(candidate, fj64_bases[fj64_hash(candidate)]);
}

/* Profiled solution finder */
static uint64_t find_solution_profiled(uint64_t n) {
    uint64_t t0, t1, t_loop;
    uint64_t t_total_start = rdtsc();

    uint64_t N = 8 * n + 3;

    /* isqrt */
    t0 = rdtsc();
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;
    t1 = rdtsc();
    cycles_isqrt += t1 - t0;

    uint64_t a = a_max;
    t_loop = rdtsc();

    while (1) {
        /* a² computation */
        t0 = rdtsc();
        uint64_t a_sq = a * a;
        t1 = rdtsc();
        cycles_a_squared += t1 - t0;

        if (a_sq <= N - 4) {
            /* candidate calculation */
            t0 = rdtsc();
            uint64_t candidate = (N - a_sq) >> 1;
            t1 = rdtsc();
            cycles_candidate_calc += t1 - t0;

            count_candidates++;

            if (is_prime_profiled(candidate)) {
                cycles_total += rdtsc() - t_total_start;
                return a;
            }
        }

        /* loop overhead */
        t0 = rdtsc();
        if (a < 3) break;
        a -= 2;
        t1 = rdtsc();
        cycles_loop_overhead += t1 - t0;
    }

    cycles_total += rdtsc() - t_total_start;
    return 0;
}

int main(void) {
    printf("Per-Line Cycle Profiling\n");
    printf("========================\n\n");

    uint64_t n_start = 1000000000000000000ULL;
    uint64_t count = 50000;

    printf("Testing %llu values starting at n = %.2e\n\n",
           (unsigned long long)count, (double)n_start);

    /* Reset counters */
    cycles_isqrt = cycles_loop_overhead = cycles_a_squared = 0;
    cycles_candidate_calc = cycles_trial_div = 0;
    cycles_montgomery_setup = cycles_montgomery_exp = cycles_montgomery_witness = 0;
    cycles_total = 0;
    count_candidates = count_trial_filtered = count_mr_tests = 0;

    for (uint64_t i = 0; i < count; i++) {
        find_solution_profiled(n_start + i);
    }

    /* Calculate percentages */
    uint64_t cycles_mr_total = cycles_montgomery_setup + cycles_montgomery_exp + cycles_montgomery_witness;
    uint64_t cycles_measured = cycles_isqrt + cycles_loop_overhead + cycles_a_squared +
                               cycles_candidate_calc + cycles_trial_div + cycles_mr_total;

    printf("=== CYCLE BREAKDOWN ===\n\n");
    printf("%-35s %12s %8s\n", "Code Section", "Cycles", "%%");
    printf("─────────────────────────────────────────────────────────\n");

    printf("%-35s %12llu %7.2f%%\n", "isqrt64(N)",
           (unsigned long long)cycles_isqrt, 100.0 * cycles_isqrt / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "a * a",
           (unsigned long long)cycles_a_squared, 100.0 * cycles_a_squared / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "(N - a_sq) >> 1",
           (unsigned long long)cycles_candidate_calc, 100.0 * cycles_candidate_calc / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "Loop control (if, a -= 2)",
           (unsigned long long)cycles_loop_overhead, 100.0 * cycles_loop_overhead / cycles_measured);

    printf("\n");
    printf("%-35s %12llu %7.2f%%\n", "Trial division (30 primes)",
           (unsigned long long)cycles_trial_div, 100.0 * cycles_trial_div / cycles_measured);

    printf("\n");
    printf("%-35s %12llu %7.2f%%\n", "Montgomery setup (n_inv, r²)",
           (unsigned long long)cycles_montgomery_setup, 100.0 * cycles_montgomery_setup / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "Montgomery exponentiation loop",
           (unsigned long long)cycles_montgomery_exp, 100.0 * cycles_montgomery_exp / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "Montgomery witness loop",
           (unsigned long long)cycles_montgomery_witness, 100.0 * cycles_montgomery_witness / cycles_measured);

    printf("─────────────────────────────────────────────────────────\n");
    printf("%-35s %12llu %7.2f%%\n", "TOTAL MEASURED",
           (unsigned long long)cycles_measured, 100.0);

    printf("\n=== GROUPED SUMMARY ===\n\n");
    uint64_t cycles_iteration = cycles_isqrt + cycles_loop_overhead + cycles_a_squared + cycles_candidate_calc;
    printf("%-35s %12llu %7.2f%%\n", "Iteration overhead (total)",
           (unsigned long long)cycles_iteration, 100.0 * cycles_iteration / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "Trial division (total)",
           (unsigned long long)cycles_trial_div, 100.0 * cycles_trial_div / cycles_measured);
    printf("%-35s %12llu %7.2f%%\n", "Miller-Rabin (total)",
           (unsigned long long)cycles_mr_total, 100.0 * cycles_mr_total / cycles_measured);

    printf("\n=== STATISTICS ===\n\n");
    printf("Candidates checked:     %llu (%.1f per n)\n",
           (unsigned long long)count_candidates, (double)count_candidates / count);
    printf("Trial-div filtered:     %llu (%.1f%% of candidates)\n",
           (unsigned long long)count_trial_filtered, 100.0 * count_trial_filtered / count_candidates);
    printf("Miller-Rabin tests:     %llu (%.1f per solution)\n",
           (unsigned long long)count_mr_tests, (double)count_mr_tests / count);

    return 0;
}
