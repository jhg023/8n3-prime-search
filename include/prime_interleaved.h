/*
 * Miller-Rabin with Interleaved FP Trial Division
 *
 * Performs trial division using floating-point reciprocals during the
 * MR exponentiation loop. FP operations run in parallel with integer
 * operations on modern CPUs (separate execution units).
 *
 * Benefits:
 * - ~4% speedup for TD-resistant candidates
 * - ~30% of composites caught by extended TD (primes 131-467)
 * - FP operations have near-zero overhead due to parallelism
 *
 * Requirements:
 * - Candidates must have already passed basic TD (primes 3-127)
 * - Candidates must be < 2^49 for FP precision (we use up to ~10^13)
 */

#ifndef PRIME_INTERLEAVED_H
#define PRIME_INTERLEAVED_H

#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include "prime.h"  /* For fj64_hash, fj64_bases, and Montgomery functions */

/* Extended primes for interleaved TD (after basic TD with 3-127) */
static const uint64_t INTERLEAVED_PRIMES[] = {
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383,
    389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467
};
#define NUM_INTERLEAVED_PRIMES 60

/* Pre-computed reciprocals with fudge factor for rounding */
static const double INTERLEAVED_RECIPROCALS[] = {
    0x1.f44659e4a44b1p-8,  /* 1/131 */
    0x1.de5d6e3f888cap-8,  /* 1/137 */
    0x1.d77b654b82e74p-8,  /* 1/139 */
    0x1.b7d6c3dda35cbp-8,  /* 1/149 */
    0x1.b2036406c8319p-8,  /* 1/151 */
    0x1.a16d3f97a4d42p-8,  /* 1/157 */
    0x1.920fb49d0e469p-8,  /* 1/163 */
    0x1.886e5f0abb28ap-8,  /* 1/167 */
    0x1.7ad2208e0ef03p-8,  /* 1/173 */
    0x1.6e1f76b433a07p-8,  /* 1/179 */
    0x1.6a13cd15374dp-8,   /* 1/181 */
    0x1.571ed3c506d7ap-8,  /* 1/191 */
    0x1.5390948f4122bp-8,  /* 1/193 */
    0x1.4cab88725b1aep-8,  /* 1/197 */
    0x1.49539e3b2d2a7p-8,  /* 1/199 */
    0x1.3698df3de0988p-8,  /* 1/211 */
    0x1.25e2270809531p-8,  /* 1/223 */
    0x1.20b470c67c319p-8,  /* 1/227 */
    0x1.1e2ef3b3fbab4p-8,  /* 1/229 */
    0x1.19453808ca4dcp-8,  /* 1/233 */
    0x1.12358e75d3273p-8,  /* 1/239 */
    0x1.0fef010fef251p-8,  /* 1/241 */
    0x1.05197f7d73644p-8,  /* 1/251 */
    0x1.fe01fe01fe4a1p-9,  /* 1/257 */
    0x1.f25f644230f36p-9,  /* 1/263 */
    0x1.e741aa5975565p-9,  /* 1/269 */
    0x1.e3a9179dc1ef4p-9,  /* 1/271 */
    0x1.d92f2231e840bp-9,  /* 1/277 */
    0x1.d272ca3fc5f9bp-9,  /* 1/281 */
    0x1.cf26e5c44c447p-9,  /* 1/283 */
    0x1.bf583ee86920cp-9,  /* 1/293 */
    0x1.aaf1d2f87f07ep-9,  /* 1/307 */
    0x1.a574107688ecbp-9,  /* 1/311 */
    0x1.a2c2a87c52121p-9,  /* 1/313 */
    0x1.9d79f176b6caep-9,  /* 1/317 */
    0x1.8bfce806303bbp-9,  /* 1/331 */
    0x1.84f00c2780a95p-9,  /* 1/337 */
    0x1.79baa6bb63e0cp-9,  /* 1/347 */
    0x1.77908119aca8ep-9,  /* 1/349 */
    0x1.734f0c542030ep-9,  /* 1/353 */
    0x1.6d1a62681cce2p-9,  /* 1/359 */
    0x1.6524f853b4f24p-9,  /* 1/367 */
    0x1.5f6643429327dp-9,  /* 1/373 */
    0x1.59d61f123d12bp-9,  /* 1/379 */
    0x1.56397ba7c5763p-9,  /* 1/383 */
    0x1.50f22e111c946p-9,  /* 1/389 */
    0x1.4a27fad7605cbp-9,  /* 1/397 */
    0x1.46dce345964e7p-9,  /* 1/401 */
    0x1.40782d10e69e7p-9,  /* 1/409 */
    0x1.38d22d3660d0fp-9,  /* 1/419 */
    0x1.3755bd1c94a6fp-9,  /* 1/421 */
    0x1.301c82ac406e1p-9,  /* 1/431 */
    0x1.2eb4ea1fed5ccp-9,  /* 1/433 */
    0x1.2a91c92f3c586p-9,  /* 1/439 */
    0x1.27dfa38a1d2cep-9,  /* 1/443 */
    0x1.23eb7971764dcp-9,  /* 1/449 */
    0x1.1ecf43c7fbccdp-9,  /* 1/457 */
    0x1.1c522fc1ce4dap-9,  /* 1/461 */
    0x1.1b17c67f2bf64p-9,  /* 1/463 */
    0x1.18ab08390305cp-9,  /* 1/467 */
};

/**
 * Check divisibility using FP reciprocal multiplication
 * Returns true if n is divisible by the prime at idx
 */
static inline bool check_fp_divisible(double n_dbl, int idx) {
    double product = n_dbl * INTERLEAVED_RECIPROCALS[idx];
    double frac = product - floor(product);
    return frac <= INTERLEAVED_RECIPROCALS[idx];
}

/**
 * Miller-Rabin witness test with interleaved FP trial division
 *
 * Performs TD checks using FP reciprocals during the exponentiation loop.
 * FP operations run in parallel with integer Montgomery operations.
 *
 * @param n         The candidate (must be odd, > 127, passed basic TD)
 * @param a         The witness base
 * @param n_inv     Pre-computed Montgomery inverse
 * @param r_sq      Pre-computed r^2 mod n
 * @param n_dbl     n as a double (for FP operations)
 * @param prime_idx Pointer to current prime index (updated on return)
 * @return          true if n passes this witness test, false if composite
 */
static inline bool mr_witness_interleaved(uint64_t n, uint64_t a,
                                          uint64_t n_inv, uint64_t r_sq,
                                          double n_dbl, int* prime_idx) {
    if (a >= n) a %= n;
    if (a == 0) return true;

    uint64_t d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;

    /* Montgomery form constants */
    uint64_t one_m = montgomery_reduce(r_sq, n, n_inv);
    uint64_t neg_one_m = n - one_m;
    uint64_t a_m = montgomery_reduce((__uint128_t)a * r_sq, n, n_inv);
    uint64_t x_m = one_m;
    uint64_t base_m = a_m;
    uint64_t exp = d;

    int idx = *prime_idx;

    /* Exponentiation loop with interleaved FP TD */
    while (exp > 0) {
        /* Check FP result from previous iteration */
        if (idx > 0 && idx <= NUM_INTERLEAVED_PRIMES) {
            if (check_fp_divisible(n_dbl, idx - 1)) {
                *prime_idx = idx;
                return false;  /* Composite: factor found */
            }
        }

        /* Montgomery exponentiation step */
        uint64_t temp = montgomery_mul(x_m, base_m, n, n_inv);
        x_m = (exp & 1) ? temp : x_m;
        base_m = montgomery_mul(base_m, base_m, n, n_inv);
        exp >>= 1;

        /* Advance to next prime for parallel FP check */
        if (idx < NUM_INTERLEAVED_PRIMES) {
            idx++;
        }
    }

    /* Final FP check */
    if (idx > 0 && idx <= NUM_INTERLEAVED_PRIMES) {
        if (check_fp_divisible(n_dbl, idx - 1)) {
            *prime_idx = idx;
            return false;
        }
    }

    /* Standard MR check */
    if (x_m == one_m || x_m == neg_one_m) {
        *prime_idx = idx;
        return true;
    }

    /* Squaring loop with continued FP TD */
    for (int i = 1; i < r; i++) {
        if (idx < NUM_INTERLEAVED_PRIMES) {
            if (check_fp_divisible(n_dbl, idx - 1)) {
                *prime_idx = idx;
                return false;
            }
            idx++;
        }

        x_m = montgomery_mul(x_m, x_m, n, n_inv);
        if (x_m == neg_one_m) {
            *prime_idx = idx;
            return true;
        }
        if (x_m == one_m) {
            *prime_idx = idx;
            return false;
        }
    }

    *prime_idx = idx;
    return false;
}

/**
 * Deterministic primality test using FJ64 with interleaved TD
 *
 * Uses 2 MR witnesses (base 2 + hash-selected) with FP trial division
 * running in parallel during the exponentiation loops.
 *
 * @param n  The candidate (must be odd, > 127, already passed basic TD)
 * @return   true if prime, false if composite
 */
static inline bool is_prime_fj64_interleaved(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    double n_dbl = (double)n;
    int prime_idx = 0;

    /* First witness: base 2 */
    if (!mr_witness_interleaved(n, 2, n_inv, r_sq, n_dbl, &prime_idx)) {
        return false;
    }

    /* Second witness: hash-selected from FJ64 table */
    uint64_t witness2 = fj64_bases[fj64_hash(n)];
    if (!mr_witness_interleaved(n, witness2, n_inv, r_sq, n_dbl, &prime_idx)) {
        return false;
    }

    return true;
}

#endif /* PRIME_INTERLEAVED_H */
