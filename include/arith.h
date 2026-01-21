/*
 * Core Arithmetic Utilities
 *
 * Modular arithmetic and integer square root functions.
 */

#ifndef ARITH_H
#define ARITH_H

#include <stdint.h>
#include <math.h>

/**
 * 64-bit modular multiplication using 128-bit intermediate
 */
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    return ((__uint128_t)a * b) % m;
}

/**
 * Modular exponentiation: base^exp mod m
 */
static inline uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1)
            result = mulmod64(result, base, mod);
        exp >>= 1;
        base = mulmod64(base, base, mod);
    }
    return result;
}

/**
 * Integer square root using Newton's method with floating-point seed
 */
static inline uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = (uint64_t)sqrtl((long double)n);
    /* Correct for floating-point errors */
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

#endif /* ARITH_H */
