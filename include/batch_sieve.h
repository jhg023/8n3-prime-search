/*
 * Batch Processing with Segmented Sieve
 *
 * For a batch of n values [n_start, n_start + batch_size), we can exploit
 * that for fixed 'a', the candidate p values form an arithmetic progression:
 *
 *   p(n) = (8n + 3 - a^2) / 2 = 4n + (3 - a^2)/2
 *
 * Instead of testing each p individually with Miller-Rabin, we can sieve
 * out composites in the progression, then only verify remaining candidates.
 *
 * This approach is better suited for exhaustive verification (checking ALL n
 * in a range) rather than early-exit per-n solution finding.
 *
 * Memory per batch (batch_size = 64K):
 *   - composite bitmap: 8 KB
 *   - solved array: 64 KB (1 byte per n)
 *   - solutions_a: 512 KB (8 bytes per n)
 *   - Total: ~584 KB per instance
 */

#ifndef BATCH_SIEVE_H
#define BATCH_SIEVE_H

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "arith.h"
#include "prime.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Default batch size (number of n values processed together) */
#define BATCH_DEFAULT_SIZE 65536

/* Sieving primes: use small primes for filtering composites in progression */
/* We use primes up to sqrt(max_p) for the batch */
#define BATCH_SIEVE_PRIME_COUNT 1000  /* Primes up to ~8000 */

/* Small primes for sieving */
static const uint32_t BATCH_SMALL_PRIMES[] = {
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151,
    157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317,
    331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503,
    509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701,
    709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811,
    821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887
};
#define BATCH_SMALL_PRIMES_COUNT 154

/* ========================================================================== */
/* Data Structures                                                            */
/* ========================================================================== */

/**
 * Batch sieve context for processing a range of n values.
 */
typedef struct {
    /* Range being processed */
    uint64_t n_start;           /* First n in batch */
    uint64_t batch_size;        /* Number of n values */

    /* Result arrays (one entry per n in batch) */
    uint8_t *solved;            /* solved[i] = 1 iff n_start + i has solution */
    uint64_t *solutions_a;      /* solutions_a[i] = the 'a' that solved it */
    uint64_t *solutions_p;      /* solutions_p[i] = the 'p' that solved it */

    /* Composite bitmap for current sieving pass */
    /* Bit i is set iff candidate for n_start + i is composite */
    uint64_t *composite;
    uint64_t composite_words;

    /* Statistics */
    uint64_t total_solved;
    uint64_t mr_tests_saved;    /* Miller-Rabin tests avoided via sieving */
    uint64_t mr_tests_done;     /* Miller-Rabin tests still needed */
} BatchSieve;

/* ========================================================================== */
/* Creation / Destruction                                                     */
/* ========================================================================== */

/**
 * Create a batch sieve context.
 * Returns NULL on allocation failure.
 */
static inline BatchSieve* batch_sieve_create(uint64_t n_start, uint64_t batch_size) {
    BatchSieve *bs = (BatchSieve*)malloc(sizeof(BatchSieve));
    if (!bs) return NULL;

    bs->n_start = n_start;
    bs->batch_size = batch_size;
    bs->total_solved = 0;
    bs->mr_tests_saved = 0;
    bs->mr_tests_done = 0;

    /* Allocate result arrays */
    bs->solved = (uint8_t*)calloc(batch_size, sizeof(uint8_t));
    bs->solutions_a = (uint64_t*)calloc(batch_size, sizeof(uint64_t));
    bs->solutions_p = (uint64_t*)calloc(batch_size, sizeof(uint64_t));

    /* Allocate composite bitmap */
    bs->composite_words = (batch_size + 63) / 64;
    bs->composite = (uint64_t*)malloc(bs->composite_words * sizeof(uint64_t));

    if (!bs->solved || !bs->solutions_a || !bs->solutions_p || !bs->composite) {
        if (bs->solved) free(bs->solved);
        if (bs->solutions_a) free(bs->solutions_a);
        if (bs->solutions_p) free(bs->solutions_p);
        if (bs->composite) free(bs->composite);
        free(bs);
        return NULL;
    }

    return bs;
}

/**
 * Destroy a batch sieve context.
 */
static inline void batch_sieve_destroy(BatchSieve *bs) {
    if (bs) {
        free(bs->solved);
        free(bs->solutions_a);
        free(bs->solutions_p);
        free(bs->composite);
        free(bs);
    }
}

/**
 * Reset the batch sieve for a new range.
 */
static inline void batch_sieve_reset(BatchSieve *bs, uint64_t n_start) {
    bs->n_start = n_start;
    bs->total_solved = 0;
    bs->mr_tests_saved = 0;
    bs->mr_tests_done = 0;

    memset(bs->solved, 0, bs->batch_size * sizeof(uint8_t));
    memset(bs->solutions_a, 0, bs->batch_size * sizeof(uint64_t));
    memset(bs->solutions_p, 0, bs->batch_size * sizeof(uint64_t));
}

/* ========================================================================== */
/* Bitmap Operations                                                          */
/* ========================================================================== */

static inline void batch_set_composite(BatchSieve *bs, uint64_t idx) {
    if (idx < bs->batch_size) {
        bs->composite[idx >> 6] |= (1ULL << (idx & 63));
    }
}

static inline bool batch_is_composite(const BatchSieve *bs, uint64_t idx) {
    if (idx >= bs->batch_size) return true;
    return (bs->composite[idx >> 6] & (1ULL << (idx & 63))) != 0;
}

static inline void batch_clear_composite_bitmap(BatchSieve *bs) {
    memset(bs->composite, 0, bs->composite_words * sizeof(uint64_t));
}

/* ========================================================================== */
/* Sieving for Fixed 'a' Value                                                */
/* ========================================================================== */

/**
 * For fixed 'a', sieve the arithmetic progression of candidate p values.
 *
 * The candidates are: p[i] = (8*(n_start + i) + 3 - a^2) / 2
 *                         = 4*(n_start + i) + c_a
 * where c_a = (3 - a^2) / 2
 *
 * These form an arithmetic progression with common difference 4.
 *
 * To sieve: for each small prime q, mark positions where p[i] % q == 0.
 * This happens every q positions (offset by first occurrence).
 */
static inline void batch_sieve_for_a(BatchSieve *bs, uint64_t a) {
    uint64_t a_sq = a * a;

    /* First candidate: p[0] = (8*n_start + 3 - a^2) / 2 */
    uint64_t N_start = 8 * bs->n_start + 3;
    if (a_sq >= N_start) {
        /* All candidates would be negative/zero - nothing to sieve */
        return;
    }
    uint64_t p_start = (N_start - a_sq) / 2;

    /* Clear composite bitmap */
    batch_clear_composite_bitmap(bs);

    /* Sieve with small primes */
    for (int i = 0; i < BATCH_SMALL_PRIMES_COUNT; i++) {
        uint32_t q = BATCH_SMALL_PRIMES[i];

        /* Find first position where p[idx] % q == 0 */
        /* p[idx] = p_start + 4*idx */
        /* We need p_start + 4*idx = 0 (mod q) */
        /* => 4*idx = -p_start (mod q) */
        /* => idx = (-p_start / 4) (mod q) */

        /* Compute p_start mod q */
        uint64_t p_mod = p_start % q;

        /* We need 4*idx = (q - p_mod) mod q */
        /* If q == 2, special case: all p values are odd (after division), skip */
        if (q == 2) continue;  /* p = (N - a^2)/2 is always odd since N is odd and a is odd */

        /* Find modular inverse of 4 mod q */
        /* For small q, we can compute 4^(-1) mod q directly */
        uint64_t inv4 = 0;
        for (uint32_t t = 0; t < q; t++) {
            if ((4 * t) % q == 1) {
                inv4 = t;
                break;
            }
        }

        /* First position: idx = inv4 * (q - p_mod) mod q */
        uint64_t first_idx = (inv4 * ((q - p_mod) % q)) % q;

        /* Mark all positions idx = first_idx + k*q as composite */
        for (uint64_t idx = first_idx; idx < bs->batch_size; idx += q) {
            /* Don't mark if this n is already solved */
            if (!bs->solved[idx]) {
                /* Check if this is actually the prime q itself */
                uint64_t p_val = p_start + 4 * idx;
                if (p_val != q) {
                    batch_set_composite(bs, idx);
                    bs->mr_tests_saved++;
                }
            }
        }
    }
}

/**
 * After sieving for 'a', check remaining candidates with Miller-Rabin.
 * Mark solutions found.
 */
static inline void batch_check_remaining(BatchSieve *bs, uint64_t a) {
    uint64_t a_sq = a * a;

    for (uint64_t idx = 0; idx < bs->batch_size; idx++) {
        /* Skip already solved */
        if (bs->solved[idx]) continue;

        /* Skip composites found by sieving */
        if (batch_is_composite(bs, idx)) continue;

        /* Compute candidate p = (8*(n_start + idx) + 3 - a^2) / 2 */
        uint64_t n = bs->n_start + idx;
        uint64_t N = 8 * n + 3;

        /* Skip if a^2 >= N (no valid p) */
        if (a_sq >= N) continue;

        uint64_t p = (N - a_sq) / 2;

        /* Skip if p < 2 */
        if (p < 2) continue;

        bs->mr_tests_done++;

        /* Test with Miller-Rabin */
        if (is_prime_fj64_fast(p)) {
            /* Solution found! */
            bs->solved[idx] = 1;
            bs->solutions_a[idx] = a;
            bs->solutions_p[idx] = p;
            bs->total_solved++;
        }
    }
}

/* ========================================================================== */
/* Full Batch Processing                                                      */
/* ========================================================================== */

/**
 * Process the entire batch, iterating through 'a' values largest-first.
 * This matches the strategy of the main search.
 */
static inline void batch_process(BatchSieve *bs) {
    /* Find a_max for the largest n in the batch */
    uint64_t n_max = bs->n_start + bs->batch_size - 1;
    uint64_t N_max = 8 * n_max + 3;
    uint64_t a_max = isqrt64(N_max);

    /* Ensure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    /* Iterate through a values, largest first */
    for (uint64_t a = a_max; a >= 1; a -= 2) {
        /* Sieve for this 'a' value */
        batch_sieve_for_a(bs, a);

        /* Check remaining candidates */
        batch_check_remaining(bs, a);

        /* Check if all solved (early termination) */
        if (bs->total_solved >= bs->batch_size) {
            break;
        }
    }
}

/**
 * Verify that all n in the batch have solutions.
 * Returns true if all solved, false if any counterexample found.
 */
static inline bool batch_verify_complete(const BatchSieve *bs) {
    for (uint64_t idx = 0; idx < bs->batch_size; idx++) {
        if (!bs->solved[idx]) {
            return false;
        }
    }
    return true;
}

/**
 * Get the number of unsolved n values (potential counterexamples).
 */
static inline uint64_t batch_count_unsolved(const BatchSieve *bs) {
    uint64_t count = 0;
    for (uint64_t idx = 0; idx < bs->batch_size; idx++) {
        if (!bs->solved[idx]) {
            count++;
        }
    }
    return count;
}

/* ========================================================================== */
/* Statistics Reporting                                                       */
/* ========================================================================== */

/**
 * Print batch statistics.
 */
static inline void batch_print_stats(const BatchSieve *bs) {
    printf("Batch Statistics:\n");
    printf("  Range: n in [%llu, %llu)\n",
           (unsigned long long)bs->n_start,
           (unsigned long long)(bs->n_start + bs->batch_size));
    printf("  Batch size: %llu\n", (unsigned long long)bs->batch_size);
    printf("  Solved: %llu / %llu\n",
           (unsigned long long)bs->total_solved,
           (unsigned long long)bs->batch_size);
    printf("  MR tests saved by sieving: %llu\n",
           (unsigned long long)bs->mr_tests_saved);
    printf("  MR tests performed: %llu\n",
           (unsigned long long)bs->mr_tests_done);

    if (bs->mr_tests_saved + bs->mr_tests_done > 0) {
        double save_rate = 100.0 * bs->mr_tests_saved /
                           (bs->mr_tests_saved + bs->mr_tests_done);
        printf("  MR test savings: %.1f%%\n", save_rate);
    }

    /* Report any unsolved (potential counterexamples) */
    uint64_t unsolved = batch_count_unsolved(bs);
    if (unsolved > 0) {
        printf("  UNSOLVED (potential counterexamples): %llu\n",
               (unsigned long long)unsolved);
        printf("  First few unsolved n values:\n");
        int shown = 0;
        for (uint64_t idx = 0; idx < bs->batch_size && shown < 5; idx++) {
            if (!bs->solved[idx]) {
                printf("    n = %llu\n", (unsigned long long)(bs->n_start + idx));
                shown++;
            }
        }
    }
}

#endif /* BATCH_SIEVE_H */
