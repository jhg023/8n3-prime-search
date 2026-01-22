/*
 * Solution Finding Utilities
 *
 * Core algorithms for finding solutions to 8n + 3 = a^2 + 2p
 *
 * Includes multiple iteration strategies:
 * - find_solution():              Large-to-small (default, fastest)
 * - find_solution_small_to_large: Original strategy
 * - find_solution_large_to_small: Alias for find_solution()
 * - find_solution_middle_out:     Start from middle, alternate outward
 * - find_solution_outside_in:     Alternate between extremes
 * - find_solution_random:         Pseudo-random order
 */

#ifndef SOLVE_H
#define SOLVE_H

#include <stdint.h>
#include <stdbool.h>
#include "arith.h"
#include "prime.h"

/* ========================================================================== */
/* Statistics Tracking (optional)                                             */
/* ========================================================================== */

/* Define SOLVE_TRACK_STATS before including this header to enable stats */
#ifdef SOLVE_TRACK_STATS
static uint64_t solve_stat_n_processed = 0;      /* Number of n values processed */
static uint64_t solve_stat_total_checks = 0;     /* Total a values checked */
static uint64_t solve_stat_candidates_32bit = 0; /* Candidates fitting in 32-bit */

static inline void solve_reset_stats(void) {
    solve_stat_n_processed = 0;
    solve_stat_total_checks = 0;
    solve_stat_candidates_32bit = 0;
}

static inline void solve_get_stats(uint64_t* n_processed, uint64_t* total_checks, uint64_t* bits32) {
    if (n_processed) *n_processed = solve_stat_n_processed;
    if (total_checks) *total_checks = solve_stat_total_checks;
    if (bits32) *bits32 = solve_stat_candidates_32bit;
}

static inline double solve_get_avg_checks(void) {
    return (solve_stat_n_processed > 0)
        ? (double)solve_stat_total_checks / solve_stat_n_processed
        : 0.0;
}

#define SOLVE_TRACK_CANDIDATE(candidate) do { \
    solve_stat_total_checks++; \
    if ((candidate) <= UINT32_MAX) solve_stat_candidates_32bit++; \
} while(0)

#define SOLVE_TRACK_SOLUTION_FOUND() do { \
    solve_stat_n_processed++; \
} while(0)
#else
#define SOLVE_TRACK_CANDIDATE(candidate) ((void)0)
#define SOLVE_TRACK_SOLUTION_FOUND() ((void)0)
#endif

/* ========================================================================== */
/* Core Trial Division Helper                                                 */
/* ========================================================================== */

/**
 * Check if candidate is filtered by trial division
 * Returns: 0 = composite (filtered), 1 = is small prime, 2 = needs Miller-Rabin
 *
 * Optimizations:
 * - Inline first 7 primes (3,5,7,11,13,17,19) which catch ~65% of composites
 * - Unroll remaining loop 4x to reduce branch overhead
 */
static inline int trial_division_check(uint64_t candidate) {
    /* Inline first 7 primes - catch ~65% of composites */
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    if (candidate % 11 == 0) return (candidate == 11) ? 1 : 0;
    if (candidate % 13 == 0) return (candidate == 13) ? 1 : 0;
    if (candidate % 17 == 0) return (candidate == 17) ? 1 : 0;
    if (candidate % 19 == 0) return (candidate == 19) ? 1 : 0;

    /* Check remaining primes 23-127 with 4x unrolled loop */
    /* Primes[7..29] = 23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127 */
    int i = 7;
    while (i + 3 < NUM_TRIAL_PRIMES) {
        if (candidate % TRIAL_PRIMES[i] == 0)
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+1] == 0)
            return (candidate == TRIAL_PRIMES[i+1]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+2] == 0)
            return (candidate == TRIAL_PRIMES[i+2]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+3] == 0)
            return (candidate == TRIAL_PRIMES[i+3]) ? 1 : 0;
        i += 4;
    }
    /* Handle remaining 1-3 primes */
    while (i < NUM_TRIAL_PRIMES) {
        if (candidate % TRIAL_PRIMES[i] == 0)
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        i++;
    }
    return 2;
}

/**
 * Test if a candidate prime is actually prime
 */
static inline bool is_candidate_prime(uint64_t candidate) {
    int td = trial_division_check(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

/* ========================================================================== */
/* Main Solution Finder (Large-to-Small - Default)                            */
/* ========================================================================== */

/**
 * Find a solution to 8n + 3 = a^2 + 2p (pre-computed N and a_max version)
 *
 * This version takes N and a_max directly, avoiding redundant isqrt64 calls
 * when the caller already has these values (e.g., in search loops where
 * a_max changes very rarely as n increases).
 *
 * Optimizations:
 * 1. Track candidate incrementally instead of computing a*a each time.
 *    When a decreases by 2: new_candidate = old_candidate + delta
 *    where delta = 2*(a-1) initially and decreases by 4 each step.
 * 2. Track delta separately to avoid multiplication in inner loop.
 *
 * Returns the largest valid a, or 0 if no solution exists (counterexample).
 * Optionally returns the corresponding prime p via out parameter.
 */
static inline uint64_t find_solution_from_N(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    /* Initialize candidate = (N - a_max²) / 2 and delta = 2*(a_max - 1) */
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    /* Iterate in reverse: largest a first means smallest candidate p first */
    while (1) {
        /* Ensure candidate p >= 2 */
        if (candidate >= 2) {
            SOLVE_TRACK_CANDIDATE(candidate);

            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                SOLVE_TRACK_SOLUTION_FOUND();
                return a;
            }
        }

        /* Decrement a by 2, checking for underflow */
        if (a < 3) break;

        /* Update candidate and delta incrementally */
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    SOLVE_TRACK_SOLUTION_FOUND();  /* Count even if no solution (counterexample) */
    return 0;  /* Counterexample! */
}

/**
 * Find a solution to 8n + 3 = a^2 + 2p
 *
 * Convenience wrapper that computes N and a_max from n.
 * For high-throughput loops, prefer find_solution_from_N() with
 * incrementally-maintained N and a_max values.
 */
static inline uint64_t find_solution(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Ensure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    return find_solution_from_N(N, a_max, p_out);
}

/* ========================================================================== */
/* Alternative Strategies (for benchmarking)                                  */
/* ========================================================================== */

/**
 * Strategy: Smallest to largest (original)
 * a = 1, 3, 5, 7, ...
 * Produces largest p values first
 */
static inline uint64_t find_solution_small_to_large(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    for (uint64_t a = 1; a <= a_max; a += 2) {
        uint64_t a_sq = a * a;

        /* Ensure candidate p >= 2 */
        if (a_sq > N - 4) break;

        uint64_t candidate = (N - a_sq) >> 1;

        SOLVE_TRACK_CANDIDATE(candidate);

        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }

    return 0;  /* Counterexample! */
}

/**
 * Strategy: Largest to smallest (alias for find_solution)
 * a = a_max, a_max-2, a_max-4, ...
 * Produces smallest p values first (small primes are more common)
 */
static inline uint64_t find_solution_large_to_small(uint64_t n, uint64_t* p_out) {
    return find_solution(n, p_out);
}

/**
 * Strategy: Middle-out
 * Start from middle a, alternate larger/smaller
 * Balances between small and large p values
 */
static inline uint64_t find_solution_middle_out(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Make sure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    /* Find actual max where p >= 2 */
    while (a_max >= 1 && a_max * a_max > N - 4) {
        a_max -= 2;
    }

    if (a_max < 1) return 0;

    /* Start from middle */
    uint64_t mid = (a_max / 2) | 1;  /* Make odd */
    uint64_t lo = mid;
    uint64_t hi = mid + 2;

    /* Test middle first */
    uint64_t candidate = (N - mid * mid) >> 1;
    SOLVE_TRACK_CANDIDATE(candidate);
    if (is_candidate_prime(candidate)) {
        if (p_out) *p_out = candidate;
        return mid;
    }

    while (lo >= 1 || hi <= a_max) {
        if (hi <= a_max) {
            candidate = (N - hi * hi) >> 1;
            SOLVE_TRACK_CANDIDATE(candidate);
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return hi;
            }
            hi += 2;
        }
        if (lo >= 1 && lo != mid) {
            candidate = (N - lo * lo) >> 1;
            SOLVE_TRACK_CANDIDATE(candidate);
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return lo;
            }
        }
        if (lo < 2) break;
        lo -= 2;
    }
    return 0;
}

/**
 * Strategy: Outside-in (alternate between smallest and largest)
 * a = 1, a_max, 3, a_max-2, 5, a_max-4, ...
 */
static inline uint64_t find_solution_outside_in(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Make sure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    /* Find actual max where p >= 2 */
    while (a_max >= 1 && a_max * a_max > N - 4) {
        a_max -= 2;
    }

    uint64_t lo = 1;
    uint64_t hi = a_max;

    while (lo <= hi) {
        /* Try small a (large p) */
        uint64_t candidate = (N - lo * lo) >> 1;
        SOLVE_TRACK_CANDIDATE(candidate);
        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return lo;
        }
        if (lo == hi) break;

        /* Try large a (small p) */
        candidate = (N - hi * hi) >> 1;
        SOLVE_TRACK_CANDIDATE(candidate);
        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return hi;
        }

        lo += 2;
        hi -= 2;
    }
    return 0;
}

/**
 * Strategy: Random order using simple LCG
 * Pseudo-random iteration to avoid systematic bias
 */
static inline uint64_t find_solution_random(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Make sure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    /* Find actual max where p >= 2 */
    while (a_max >= 1 && a_max * a_max > N - 4) {
        a_max -= 2;
    }

    /* Number of odd values from 1 to a_max */
    uint64_t count = (a_max + 1) / 2;
    if (count == 0) return 0;

    /* Use n as seed for reproducibility */
    uint64_t state = n * 6364136223846793005ULL + 1442695040888963407ULL;

    for (uint64_t i = 0; i < count; i++) {
        /* Generate pseudo-random index */
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        uint64_t idx = (state >> 33) % count;
        uint64_t a = 2 * idx + 1;

        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;

        uint64_t candidate = (N - a_sq) >> 1;
        SOLVE_TRACK_CANDIDATE(candidate);
        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }

    /* Fallback: linear scan to ensure we don't miss any */
    for (uint64_t a = 1; a <= a_max; a += 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) break;

        uint64_t candidate = (N - a_sq) >> 1;
        SOLVE_TRACK_CANDIDATE(candidate);
        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }
    return 0;
}

#endif /* SOLVE_H */
