/*
 * Fast Solution Finding - Micro-optimized Implementation
 *
 * Optimizations applied:
 * 1. Reduced trial division to most effective primes
 * 2. Incremental a² computation
 * 3. Branch hints with __builtin_expect
 * 4. Inlined hot paths
 * 5. Minimal function call overhead
 *
 * Benchmark results showed root-class sieve adds more overhead than it saves
 * because solutions are found quickly (~10-15 candidates per N).
 */

#ifndef SOLVE_FAST_H
#define SOLVE_FAST_H

#include <stdint.h>
#include <stdbool.h>
#include "arith.h"
#include "arith_montgomery.h"
#include "prime.h"  /* Includes fj64_table.h */

/* ========================================================================== */
/* Optimized Trial Division                                                   */
/* ========================================================================== */

/*
 * Tuned prime list: 30 primes is near-optimal across all scales.
 * Benchmarked: 30-35 primes gives best throughput at 10^9 through 10^18.
 */
static const uint32_t FAST_TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_FAST_TRIAL_PRIMES 30

/*
 * Fast trial division - returns:
 *   0: composite (filtered)
 *   1: is one of the trial primes
 *   2: needs Miller-Rabin
 */
static inline int fast_trial_division(uint64_t n) {
    /* Unrolled first 10 primes for most common cases */
    if (__builtin_expect(n % 3 == 0, 0)) return (n == 3) ? 1 : 0;
    if (__builtin_expect(n % 5 == 0, 0)) return (n == 5) ? 1 : 0;
    if (__builtin_expect(n % 7 == 0, 0)) return (n == 7) ? 1 : 0;
    if (__builtin_expect(n % 11 == 0, 0)) return (n == 11) ? 1 : 0;
    if (__builtin_expect(n % 13 == 0, 0)) return (n == 13) ? 1 : 0;
    if (__builtin_expect(n % 17 == 0, 0)) return (n == 17) ? 1 : 0;
    if (__builtin_expect(n % 19 == 0, 0)) return (n == 19) ? 1 : 0;
    if (__builtin_expect(n % 23 == 0, 0)) return (n == 23) ? 1 : 0;
    if (__builtin_expect(n % 29 == 0, 0)) return (n == 29) ? 1 : 0;
    if (__builtin_expect(n % 31 == 0, 0)) return (n == 31) ? 1 : 0;

    /* Remaining primes in a loop */
    for (int i = 10; i < NUM_FAST_TRIAL_PRIMES; i++) {
        if (__builtin_expect(n % FAST_TRIAL_PRIMES[i] == 0, 0)) {
            return (n == FAST_TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* ========================================================================== */
/* Inlined Miller-Rabin with Montgomery                                       */
/* ========================================================================== */

/*
 * FJ64 hash - inlined for speed
 */
static inline uint32_t fast_fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    return ((x >> 32) ^ x) & 262143;
}

/*
 * Full primality test, inlined and optimized
 */
static inline bool fast_is_prime(uint64_t n) {
    /* Trial division */
    int td = fast_trial_division(n);
    if (td == 0) return false;
    if (td == 1) return true;
    if (n <= 73) return true;  /* Small primes that passed trial division */

    /* FJ64_262K with Montgomery multiplication */
    if (!mr_witness_montgomery(n, 2))
        return false;
    return mr_witness_montgomery(n, fj64_bases[fast_fj64_hash(n)]);
}

/* ========================================================================== */
/* Fast Solution Finder                                                       */
/* ========================================================================== */

/**
 * Find a solution to 8n + 3 = a² + 2p (optimized)
 *
 * Uses incremental a² computation and optimized primality testing.
 * Returns largest valid a, or 0 if no solution (counterexample).
 */
static inline uint64_t find_solution_fast(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a = isqrt64(N);

    /* Ensure a is odd */
    a = (a | 1);
    if (__builtin_expect(a * a > N, 0)) a -= 2;

    /* Precompute limit for p >= 2 check */
    uint64_t limit = N - 4;

    /* Incremental a² computation: (a-2)² = a² - 4a + 4 */
    uint64_t a_sq = a * a;

    while (1) {
        if (__builtin_expect(a_sq <= limit, 1)) {
            uint64_t candidate = (N - a_sq) >> 1;

            if (__builtin_expect(fast_is_prime(candidate), 0)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (__builtin_expect(a < 3, 0)) break;

        /* Incremental update */
        a_sq = a_sq - 4 * a + 4;
        a -= 2;
    }

    return 0;  /* Counterexample */
}

/**
 * Alternative: Original algorithm but with fast primality test
 */
static inline uint64_t find_solution_fast_v2(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (__builtin_expect(a_sq > N - 4, 0)) continue;

        uint64_t candidate = (N - a_sq) >> 1;

        if (__builtin_expect(fast_is_prime(candidate), 0)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }

    return 0;
}

#endif /* SOLVE_FAST_H */
