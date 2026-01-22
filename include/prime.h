/*
 * Primality Testing Utilities
 *
 * FJ64_262K primality test and Miller-Rabin implementation.
 * Uses the Forisek-Jancina hash table for optimal witness selection.
 *
 * Optimized with Montgomery multiplication for 3x faster modular arithmetic
 * when n < 2^63 (falls back to __uint128_t for larger n).
 */

#ifndef PRIME_H
#define PRIME_H

#include <stdint.h>
#include <stdbool.h>
#include "arith.h"
#include "arith_montgomery.h"
#include "fj64_table.h"

/* ========================================================================== */
/* Trial Division Configuration                                               */
/* ========================================================================== */

/*
 * Trial division primes - 30 odd primes from 3 to 127
 * Filters ~80% of composites before Miller-Rabin
 *
 * With Montgomery multiplication (3x faster MR), fewer trial primes is optimal:
 * - 30 primes: ~80% filter rate, minimal overhead
 * - 120 primes: ~85% filter rate, but 90 extra modulo ops
 * Benchmarked: 30 primes is 7-11% faster than 120 at large n
 */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
static const int NUM_TRIAL_PRIMES = 30;

/* ========================================================================== */
/* FJ64_262K Primality Test                                                   */
/* ========================================================================== */

/**
 * FJ64 hash function - maps n to a bucket in [0, 262143]
 */
static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;  /* & (2^18 - 1) */
}

/**
 * Single Miller-Rabin witness test (standard, for backward compatibility)
 * Returns true if n is a strong probable prime to base a
 */
static inline bool mr_witness(uint64_t n, uint64_t a) {
    if (a >= n) a %= n;
    if (a == 0) return true;

    uint64_t d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;

    uint64_t x = powmod64(a, d, n);

    if (x == 1 || x == n - 1)
        return true;

    for (int i = 1; i < r; i++) {
        x = mulmod64(x, x, n);
        if (x == n - 1)
            return true;
        if (x == 1)
            return false;
    }
    return false;
}

/**
 * FJ64_262K primality test using Montgomery multiplication
 * Exactly 2 Miller-Rabin tests with optimized modular arithmetic
 * Assumes: n > 127, n is odd, n passed trial division
 *
 * Optimization: Pre-compute Montgomery constants once and reuse for both witnesses
 */
static inline bool is_prime_fj64_fast(uint64_t n) {
    /* Pre-compute Montgomery constants once for both witness tests */
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);

    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

/**
 * FJ64_262K primality test (standard version, no Montgomery)
 * For benchmarking comparison
 */
static inline bool is_prime_fj64_standard(uint64_t n) {
    if (!mr_witness(n, 2))
        return false;
    return mr_witness(n, fj64_bases[fj64_hash(n)]);
}

/**
 * Full primality test for standalone use
 */
static inline bool is_prime_64(uint64_t n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if ((n & 1) == 0) return false;
    if (n % 3 == 0) return n == 3;
    if (n % 5 == 0) return n == 5;
    if (n % 7 == 0) return n == 7;
    if (n < 121) return true;

    return is_prime_fj64_fast(n);
}

#endif /* PRIME_H */
