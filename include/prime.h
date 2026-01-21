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
 * Trial division primes - 120 odd primes from 3 to 661
 * Filters ~82% of composites before Miller-Rabin
 * Benchmarked at n ~ 2.3e18: +12.4% throughput vs 30 primes
 */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
    359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
    479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661
};
static const int NUM_TRIAL_PRIMES = 120;

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
 */
static inline bool is_prime_fj64_fast(uint64_t n) {
    if (!mr_witness_montgomery(n, 2))
        return false;
    return mr_witness_montgomery(n, fj64_bases[fj64_hash(n)]);
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
