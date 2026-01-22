/*
 * Final Optimization Benchmark
 *
 * Test the best combined optimizations vs baseline:
 * - Inline first 3 primes in trial division
 * - Defer FJ64 hash until after trial division passes
 *
 * Also verify correctness.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/fj64_table.h"

/* Trial primes */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_TRIAL_PRIMES 30

#define ITERATIONS 10000000ULL

/* ========================================================================== */
/* Current Implementation (baseline)                                          */
/* ========================================================================== */

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline int trial_division_current(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_prime_fj64_current(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

static inline bool is_candidate_prime_current(uint64_t candidate) {
    int td = trial_division_current(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_current(candidate);
}

static inline uint64_t find_solution_current(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_current(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Optimized Implementation                                                   */
/* ========================================================================== */

/**
 * Optimized trial division with inline first 3 primes
 * 3, 5, 7 catch ~50% of composites - inlining avoids loop overhead for them
 */
static inline int trial_division_optimized(uint64_t candidate) {
    /* Inline the first 3 primes (catch ~50% of composites) */
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;

    /* Check remaining primes 11-127 via loop */
    for (int i = 3; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/**
 * Optimized primality test - defer hash computation
 *
 * Key insight: ~80% of candidates are filtered by trial division.
 * The FJ64 hash is 3 multiplications - defer it until we know we need MR.
 */
static inline bool is_candidate_prime_optimized(uint64_t candidate) {
    int td = trial_division_optimized(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;

    /* Only compute Montgomery constants and hash when needed */
    uint64_t n_inv = montgomery_inverse(candidate);
    uint64_t r_sq = montgomery_r_squared(candidate);

    if (!mr_witness_montgomery_cached(candidate, 2, n_inv, r_sq))
        return false;

    /* Compute hash only for probable primes (passed base-2 test) */
    return mr_witness_montgomery_cached(candidate, fj64_bases[fj64_hash(candidate)], n_inv, r_sq);
}

static inline uint64_t find_solution_optimized(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_optimized(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }
    return 0;
}

/* ========================================================================== */
/* Verification                                                               */
/* ========================================================================== */

bool verify_correctness(uint64_t n_start, uint64_t count) {
    printf("Verifying correctness (%llu values from n=%llu)...\n",
           (unsigned long long)count, (unsigned long long)n_start);

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t mismatches = 0;

    for (uint64_t i = 0; i < count; i++) {
        uint64_t p1, p2;
        uint64_t a1 = find_solution_current(N, a_max, &p1);
        uint64_t a2 = find_solution_optimized(N, a_max, &p2);

        if (a1 != a2 || (a1 > 0 && p1 != p2)) {
            mismatches++;
            if (mismatches <= 5) {
                printf("MISMATCH at N=%llu: current=(%llu,%llu), optimized=(%llu,%llu)\n",
                       (unsigned long long)N,
                       (unsigned long long)a1, (unsigned long long)p1,
                       (unsigned long long)a2, (unsigned long long)p2);
            }
        }

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    if (mismatches == 0) {
        printf("All %llu values verified OK\n", (unsigned long long)count);
        return true;
    } else {
        printf("FAILED: %llu mismatches\n", (unsigned long long)mismatches);
        return false;
    }
}

/* ========================================================================== */
/* Benchmarks                                                                 */
/* ========================================================================== */

void benchmark(uint64_t n_start, const char* label,
               uint64_t (*find_fn)(uint64_t, uint64_t, uint64_t*)) {
    clock_t start = clock();
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        sum += find_fn(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("%-20s: %.0f n/sec (%.3f sec)\n", label, rate, elapsed);
}

int main(void) {
    printf("Final Optimization Benchmark\n\n");

    /* Verify correctness first */
    if (!verify_correctness(1000000000000ULL, 100000)) {
        return 1;
    }
    if (!verify_correctness(2000000000000000000ULL, 100000)) {
        return 1;
    }
    printf("\n");

    printf("=== Performance Comparison (%llu iterations) ===\n\n", ITERATIONS);

    printf("At n = 10^12:\n");
    benchmark(1000000000000ULL, "Current", find_solution_current);
    benchmark(1000000000000ULL, "Optimized", find_solution_optimized);

    printf("\nAt n = 10^15:\n");
    benchmark(1000000000000000ULL, "Current", find_solution_current);
    benchmark(1000000000000000ULL, "Optimized", find_solution_optimized);

    printf("\nAt n = 2e18:\n");
    benchmark(2000000000000000000ULL, "Current", find_solution_current);
    benchmark(2000000000000000000ULL, "Optimized", find_solution_optimized);

    return 0;
}
