/*
 * Prime Sieve with O(1) Lookup
 *
 * Eratosthenes sieve stored as a bitmap for memory-efficient,
 * constant-time primality testing of small numbers.
 *
 * Usage:
 *   PrimeSieve* sieve = sieve_create(100000000);  // Sieve up to 10^8
 *   if (sieve_is_prime(sieve, n)) { ... }
 *   sieve_destroy(sieve);
 *
 * Memory footprint:
 *   - 10^7: ~1.25 MB
 *   - 10^8: ~12.5 MB
 *   - 10^9: ~125 MB
 */

#ifndef PRIME_SIEVE_H
#define PRIME_SIEVE_H

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ========================================================================== */
/* Data Structures                                                            */
/* ========================================================================== */

typedef struct {
    uint64_t *bitmap;       /* Bit k set iff k is prime (odd k only after 2) */
    uint64_t threshold;     /* Maximum value in sieve */
    uint64_t num_words;     /* Number of 64-bit words in bitmap */
    uint64_t prime_count;   /* Count of primes <= threshold */
} PrimeSieve;

/* ========================================================================== */
/* Bit Manipulation Helpers                                                   */
/* ========================================================================== */

/*
 * Bitmap layout for odd numbers only (halves memory usage):
 * - Bit 0 represents 3
 * - Bit 1 represents 5
 * - Bit k represents (2k + 3)
 * - Number n (odd, >= 3) maps to bit (n - 3) / 2
 *
 * A set bit means the number is PRIME.
 */

static inline uint64_t sieve_bit_index(uint64_t n) {
    return (n - 3) >> 1;  /* (n - 3) / 2 */
}

static inline uint64_t sieve_word_index(uint64_t bit_idx) {
    return bit_idx >> 6;  /* bit_idx / 64 */
}

static inline uint64_t sieve_bit_mask(uint64_t bit_idx) {
    return 1ULL << (bit_idx & 63);  /* bit_idx % 64 */
}

static inline bool sieve_get_bit(const uint64_t *bitmap, uint64_t bit_idx) {
    return (bitmap[sieve_word_index(bit_idx)] & sieve_bit_mask(bit_idx)) != 0;
}

static inline void sieve_clear_bit(uint64_t *bitmap, uint64_t bit_idx) {
    bitmap[sieve_word_index(bit_idx)] &= ~sieve_bit_mask(bit_idx);
}

/* ========================================================================== */
/* Sieve Creation                                                             */
/* ========================================================================== */

/**
 * Create a prime sieve for all numbers up to threshold.
 * Returns NULL on allocation failure.
 *
 * Time complexity: O(n log log n) for sieve generation
 * Space complexity: O(n/16) bytes (bit per odd number)
 */
static inline PrimeSieve* sieve_create(uint64_t threshold) {
    if (threshold < 2) {
        threshold = 2;
    }

    PrimeSieve *sieve = (PrimeSieve*)malloc(sizeof(PrimeSieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;

    /* Calculate bitmap size for odd numbers >= 3 */
    uint64_t max_bit_idx = (threshold >= 3) ? sieve_bit_index(threshold | 1) : 0;
    sieve->num_words = (max_bit_idx >> 6) + 1;

    sieve->bitmap = (uint64_t*)malloc(sieve->num_words * sizeof(uint64_t));
    if (!sieve->bitmap) {
        free(sieve);
        return NULL;
    }

    /* Initialize all bits to 1 (assume prime), then sieve out composites */
    memset(sieve->bitmap, 0xFF, sieve->num_words * sizeof(uint64_t));

    /* Sieve of Eratosthenes for odd numbers */
    /* For each prime p >= 3, mark all odd multiples of p as composite */
    uint64_t sqrt_thresh = 1;
    while ((sqrt_thresh + 1) * (sqrt_thresh + 1) <= threshold) {
        sqrt_thresh++;
    }

    for (uint64_t p = 3; p <= sqrt_thresh; p += 2) {
        /* Check if p is still marked prime */
        if (!sieve_get_bit(sieve->bitmap, sieve_bit_index(p))) {
            continue;
        }

        /* Mark odd multiples of p starting from p*p */
        /* p*p is the first unmarked multiple (smaller multiples already marked) */
        /* We skip even multiples, so we step by 2*p to hit only odd multiples */
        for (uint64_t m = p * p; m <= threshold; m += 2 * p) {
            sieve_clear_bit(sieve->bitmap, sieve_bit_index(m));
        }
    }

    /* Count primes for statistics */
    sieve->prime_count = (threshold >= 2) ? 1 : 0;  /* Count 2 if in range */
    for (uint64_t n = 3; n <= threshold; n += 2) {
        if (sieve_get_bit(sieve->bitmap, sieve_bit_index(n))) {
            sieve->prime_count++;
        }
    }

    return sieve;
}

/* ========================================================================== */
/* Sieve Destruction                                                          */
/* ========================================================================== */

/**
 * Free all memory associated with the sieve.
 */
static inline void sieve_destroy(PrimeSieve *sieve) {
    if (sieve) {
        free(sieve->bitmap);
        free(sieve);
    }
}

/* ========================================================================== */
/* Primality Lookup                                                           */
/* ========================================================================== */

/**
 * Check if n is prime using O(1) bitmap lookup.
 * Requires: n <= sieve->threshold
 *
 * Returns false for n > threshold (caller should use fallback test).
 */
static inline bool sieve_is_prime(const PrimeSieve *sieve, uint64_t n) {
    if (n > sieve->threshold) {
        return false;  /* Cannot determine - caller should use fallback */
    }
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;  /* Even and > 2 */

    return sieve_get_bit(sieve->bitmap, sieve_bit_index(n));
}

/**
 * Check if n is within the sieve's range.
 */
static inline bool sieve_in_range(const PrimeSieve *sieve, uint64_t n) {
    return n <= sieve->threshold;
}

/* ========================================================================== */
/* Statistics                                                                 */
/* ========================================================================== */

/**
 * Get the number of primes in the sieve (pi(threshold)).
 */
static inline uint64_t sieve_prime_count(const PrimeSieve *sieve) {
    return sieve->prime_count;
}

/**
 * Get the memory usage of the sieve in bytes.
 */
static inline uint64_t sieve_memory_bytes(const PrimeSieve *sieve) {
    return sizeof(PrimeSieve) + sieve->num_words * sizeof(uint64_t);
}

/**
 * Get human-readable memory usage description.
 * Uses static buffer - not thread-safe.
 */
static inline const char* sieve_memory_str(const PrimeSieve *sieve) {
    static char buf[32];
    uint64_t bytes = sieve_memory_bytes(sieve);

    if (bytes < 1024) {
        snprintf(buf, sizeof(buf), "%llu B", (unsigned long long)bytes);
    } else if (bytes < 1024 * 1024) {
        snprintf(buf, sizeof(buf), "%.1f KB", bytes / 1024.0);
    } else if (bytes < 1024 * 1024 * 1024) {
        snprintf(buf, sizeof(buf), "%.1f MB", bytes / (1024.0 * 1024.0));
    } else {
        snprintf(buf, sizeof(buf), "%.2f GB", bytes / (1024.0 * 1024.0 * 1024.0));
    }
    return buf;
}

#endif /* PRIME_SIEVE_H */
