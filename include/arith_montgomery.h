/*
 * Montgomery Arithmetic for 64-bit Modular Operations
 *
 * Montgomery multiplication replaces division with multiplication and bitwise
 * operations, which is faster for repeated modular multiplications with the
 * same modulus (like in Miller-Rabin exponentiation).
 *
 * IMPORTANT: Standard Montgomery with r=2^64 can overflow for n >= 2^63.
 * For n close to 2^64, we fall back to __uint128_t division.
 *
 * Reference: https://cp-algorithms.com/algebra/montgomery_multiplication.html
 */

#ifndef ARITH_MONTGOMERY_H
#define ARITH_MONTGOMERY_H

#include <stdint.h>
#include <stdbool.h>

/* Threshold: Montgomery is safe for n < 2^63 */
#define MONTGOMERY_SAFE_THRESHOLD (1ULL << 63)

/* ========================================================================== */
/* Standard __uint128_t operations (fallback for large n)                     */
/* ========================================================================== */

static inline uint64_t mulmod64_std(uint64_t a, uint64_t b, uint64_t m) {
    return ((__uint128_t)a * b) % m;
}

static inline uint64_t powmod64_std(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1)
            result = mulmod64_std(result, base, mod);
        exp >>= 1;
        base = mulmod64_std(base, base, mod);
    }
    return result;
}

/* ========================================================================== */
/* Montgomery Arithmetic (fast path for n < 2^63)                             */
/* ========================================================================== */

/**
 * Compute n' = -n^(-1) mod 2^64 using Newton's method
 * Only valid for odd n (which is true for all our prime candidates)
 */
static inline uint64_t montgomery_inverse(uint64_t n) {
    uint64_t x = n;
    x *= 2 - n * x;  /* 4 bits */
    x *= 2 - n * x;  /* 8 bits */
    x *= 2 - n * x;  /* 16 bits */
    x *= 2 - n * x;  /* 32 bits */
    x *= 2 - n * x;  /* 64 bits */
    return -x;
}

/**
 * Montgomery reduction: compute t * r^(-1) mod n
 * ONLY safe for n < 2^63
 */
static inline uint64_t montgomery_reduce(__uint128_t t, uint64_t n, uint64_t n_inv) {
    uint64_t m = (uint64_t)t * n_inv;
    __uint128_t mn = (__uint128_t)m * n;
    uint64_t u = (t + mn) >> 64;
    return (u >= n) ? (u - n) : u;
}

/**
 * Montgomery multiplication: a * b * r^(-1) mod n
 */
static inline uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t n, uint64_t n_inv) {
    return montgomery_reduce((__uint128_t)a * b, n, n_inv);
}

/**
 * Compute r^2 mod n = 2^128 mod n
 */
static inline uint64_t montgomery_r_squared(uint64_t n) {
    uint64_t r = (((__uint128_t)1 << 64) % n);
    return ((__uint128_t)r * r) % n;
}

/* ========================================================================== */
/* Hybrid Miller-Rabin: Montgomery for small n, fallback for large n          */
/* ========================================================================== */

/**
 * Miller-Rabin witness test using Montgomery (n < 2^63)
 */
static inline bool mr_witness_montgomery_safe(uint64_t n, uint64_t a, uint64_t n_inv, uint64_t r_sq) {
    if (a >= n) a %= n;
    if (a == 0) return true;

    uint64_t d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;

    /* Montgomery form: one_m = r mod n, neg_one_m = (n-1)*r mod n */
    uint64_t one_m = montgomery_reduce(r_sq, n, n_inv);
    uint64_t neg_one_m = n - one_m;

    /* Compute a^d mod n using branchless exponentiation */
    uint64_t a_m = montgomery_reduce((__uint128_t)a * r_sq, n, n_inv);
    uint64_t x_m = one_m;
    uint64_t base_m = a_m;
    uint64_t exp = d;

    /*
     * Branchless exponentiation: always compute the multiply, use conditional
     * select. This avoids branch mispredictions which are costly when the
     * exponent bits are unpredictable.
     */
    while (exp > 0) {
        uint64_t temp = montgomery_mul(x_m, base_m, n, n_inv);
        x_m = (exp & 1) ? temp : x_m;
        base_m = montgomery_mul(base_m, base_m, n, n_inv);
        exp >>= 1;
    }

    if (x_m == one_m || x_m == neg_one_m)
        return true;

    for (int i = 1; i < r; i++) {
        x_m = montgomery_mul(x_m, x_m, n, n_inv);
        if (x_m == neg_one_m)
            return true;
        if (x_m == one_m)
            return false;
    }
    return false;
}

/**
 * Miller-Rabin witness test using standard arithmetic (fallback)
 */
static inline bool mr_witness_std(uint64_t n, uint64_t a) {
    if (a >= n) a %= n;
    if (a == 0) return true;

    uint64_t d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;

    uint64_t x = powmod64_std(a, d, n);

    if (x == 1 || x == n - 1)
        return true;

    for (int i = 1; i < r; i++) {
        x = mulmod64_std(x, x, n);
        if (x == n - 1)
            return true;
        if (x == 1)
            return false;
    }
    return false;
}

/**
 * Hybrid Miller-Rabin witness test
 * Uses Montgomery for n < 2^63, fallback otherwise
 */
static inline bool mr_witness_montgomery(uint64_t n, uint64_t a) {
    if (n < MONTGOMERY_SAFE_THRESHOLD) {
        uint64_t n_inv = montgomery_inverse(n);
        uint64_t r_sq = montgomery_r_squared(n);
        return mr_witness_montgomery_safe(n, a, n_inv, r_sq);
    } else {
        return mr_witness_std(n, a);
    }
}

/**
 * Miller-Rabin witness test with pre-computed Montgomery constants
 * Use when calling multiple witnesses on the same n
 */
static inline bool mr_witness_montgomery_cached(uint64_t n, uint64_t a,
                                                 uint64_t n_inv, uint64_t r_sq) {
    if (n < MONTGOMERY_SAFE_THRESHOLD) {
        return mr_witness_montgomery_safe(n, a, n_inv, r_sq);
    } else {
        return mr_witness_std(n, a);
    }
}

#endif /* ARITH_MONTGOMERY_H */
