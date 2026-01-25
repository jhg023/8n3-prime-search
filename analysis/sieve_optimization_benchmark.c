/*
 * Sieve Optimization Benchmark
 *
 * Tests various optimization techniques from primesieve:
 * 1. Baseline (current implementation)
 * 2. Segmented sieve (L1 cache-friendly)
 * 3. Pre-sieve pattern (multiples of 3,5,7 pre-removed)
 * 4. Mod 30 wheel (8 bits per 30 numbers)
 * 5. OpenMP parallelization
 * 6. Combined optimizations
 *
 * Usage: ./sieve_benchmark [threshold]
 * Default threshold: 10^9
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ========================================================================== */
/* Timing utilities                                                           */
/* ========================================================================== */

static double get_time_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
}

/* ========================================================================== */
/* 1. BASELINE: Current implementation from prime_sieve.h                     */
/* ========================================================================== */

typedef struct {
    uint64_t *bitmap;
    uint64_t threshold;
    uint64_t num_words;
    uint64_t prime_count;
} BaselineSieve;

static inline uint64_t baseline_bit_index(uint64_t n) {
    return (n - 3) >> 1;
}

static inline bool baseline_get_bit(const uint64_t *bitmap, uint64_t bit_idx) {
    return (bitmap[bit_idx >> 6] & (1ULL << (bit_idx & 63))) != 0;
}

static inline void baseline_clear_bit(uint64_t *bitmap, uint64_t bit_idx) {
    bitmap[bit_idx >> 6] &= ~(1ULL << (bit_idx & 63));
}

static BaselineSieve* baseline_create(uint64_t threshold) {
    BaselineSieve *sieve = malloc(sizeof(BaselineSieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;
    uint64_t max_bit_idx = baseline_bit_index(threshold | 1);
    sieve->num_words = (max_bit_idx >> 6) + 1;

    sieve->bitmap = malloc(sieve->num_words * sizeof(uint64_t));
    if (!sieve->bitmap) { free(sieve); return NULL; }

    memset(sieve->bitmap, 0xFF, sieve->num_words * sizeof(uint64_t));

    uint64_t sqrt_thresh = (uint64_t)sqrt((double)threshold) + 1;

    for (uint64_t p = 3; p <= sqrt_thresh; p += 2) {
        if (!baseline_get_bit(sieve->bitmap, baseline_bit_index(p))) continue;
        for (uint64_t m = p * p; m <= threshold; m += 2 * p) {
            baseline_clear_bit(sieve->bitmap, baseline_bit_index(m));
        }
    }

    return sieve;
}

static uint64_t baseline_count(BaselineSieve *sieve) {
    uint64_t count = (sieve->threshold >= 2) ? 1 : 0;
    for (uint64_t n = 3; n <= sieve->threshold; n += 2) {
        if (baseline_get_bit(sieve->bitmap, baseline_bit_index(n))) count++;
    }
    return count;
}

static void baseline_destroy(BaselineSieve *sieve) {
    if (sieve) { free(sieve->bitmap); free(sieve); }
}

/* ========================================================================== */
/* 2. SEGMENTED: L1 cache-friendly segmented sieve                            */
/* ========================================================================== */

#define SEGMENT_SIZE (32 * 1024)  /* 32KB L1 cache */
#define SEGMENT_BITS (SEGMENT_SIZE * 8)
#define SEGMENT_RANGE (SEGMENT_BITS * 2)  /* Each bit covers 2 numbers (odds only) */

typedef struct {
    uint64_t *bitmap;
    uint64_t threshold;
    uint64_t num_words;
} SegmentedSieve;

static SegmentedSieve* segmented_create(uint64_t threshold) {
    SegmentedSieve *sieve = malloc(sizeof(SegmentedSieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;
    uint64_t max_bit_idx = (threshold - 3) / 2;
    sieve->num_words = (max_bit_idx / 64) + 1;

    sieve->bitmap = malloc(sieve->num_words * sizeof(uint64_t));
    if (!sieve->bitmap) { free(sieve); return NULL; }
    memset(sieve->bitmap, 0xFF, sieve->num_words * sizeof(uint64_t));

    uint64_t sqrt_thresh = (uint64_t)sqrt((double)threshold) + 1;

    /* First, sieve small primes up to sqrt(threshold) */
    uint64_t small_limit = sqrt_thresh;
    uint64_t small_bits = (small_limit - 3) / 2 + 1;
    uint64_t small_words = (small_bits / 64) + 1;
    uint64_t *small_sieve = malloc(small_words * sizeof(uint64_t));
    memset(small_sieve, 0xFF, small_words * sizeof(uint64_t));

    uint64_t sqrt_small = (uint64_t)sqrt((double)small_limit) + 1;
    for (uint64_t p = 3; p <= sqrt_small; p += 2) {
        uint64_t idx = (p - 3) / 2;
        if (!((small_sieve[idx >> 6] >> (idx & 63)) & 1)) continue;
        for (uint64_t m = p * p; m <= small_limit; m += 2 * p) {
            uint64_t midx = (m - 3) / 2;
            small_sieve[midx >> 6] &= ~(1ULL << (midx & 63));
        }
    }

    /* Collect sieving primes */
    uint64_t *sieving_primes = malloc((small_limit / 2) * sizeof(uint64_t));
    uint64_t num_sieving_primes = 0;
    for (uint64_t p = 3; p <= small_limit; p += 2) {
        uint64_t idx = (p - 3) / 2;
        if ((small_sieve[idx >> 6] >> (idx & 63)) & 1) {
            sieving_primes[num_sieving_primes++] = p;
        }
    }
    free(small_sieve);

    /* Sieve in segments */
    uint8_t *segment = malloc(SEGMENT_SIZE);

    for (uint64_t seg_start = 3; seg_start <= threshold; seg_start += SEGMENT_RANGE) {
        uint64_t seg_end = seg_start + SEGMENT_RANGE - 1;
        if (seg_end > threshold) seg_end = threshold;

        memset(segment, 0xFF, SEGMENT_SIZE);

        for (uint64_t i = 0; i < num_sieving_primes; i++) {
            uint64_t p = sieving_primes[i];
            uint64_t start = p * p;
            if (start > seg_end) break;

            if (start < seg_start) {
                start = seg_start + (p - (seg_start % p)) % p;
                if ((start & 1) == 0) start += p;
            }

            for (uint64_t m = start; m <= seg_end; m += 2 * p) {
                uint64_t local_idx = (m - seg_start) / 2;
                segment[local_idx >> 3] &= ~(1 << (local_idx & 7));
            }
        }

        /* Copy segment to main bitmap */
        for (uint64_t n = seg_start; n <= seg_end; n += 2) {
            uint64_t local_idx = (n - seg_start) / 2;
            bool is_prime = (segment[local_idx >> 3] >> (local_idx & 7)) & 1;
            uint64_t global_idx = (n - 3) / 2;
            if (!is_prime) {
                sieve->bitmap[global_idx >> 6] &= ~(1ULL << (global_idx & 63));
            }
        }
    }

    free(segment);
    free(sieving_primes);
    return sieve;
}

static uint64_t segmented_count(SegmentedSieve *sieve) {
    uint64_t count = (sieve->threshold >= 2) ? 1 : 0;
    for (uint64_t n = 3; n <= sieve->threshold; n += 2) {
        uint64_t idx = (n - 3) / 2;
        if ((sieve->bitmap[idx >> 6] >> (idx & 63)) & 1) count++;
    }
    return count;
}

static void segmented_destroy(SegmentedSieve *sieve) {
    if (sieve) { free(sieve->bitmap); free(sieve); }
}

/* ========================================================================== */
/* 3. PRE-SIEVE: Use precomputed pattern for 3,5,7                            */
/* ========================================================================== */

/* Pattern for odd numbers with multiples of 3,5,7 removed */
/* LCM(3,5,7) = 105, but we only store odd numbers, so 105/2 ~= 53 bytes */
/* Actually we need the pattern to repeat correctly for odd indices */
/* For simplicity, use LCM approach: pattern repeats every 105 odd numbers */

#define PRESIEVE_PRIMES_PRODUCT 105  /* 3 * 5 * 7 */

static uint8_t presieve_pattern[PRESIEVE_PRIMES_PRODUCT];
static bool presieve_initialized = false;

static void init_presieve_pattern(void) {
    if (presieve_initialized) return;

    /* Initialize all as prime */
    memset(presieve_pattern, 0xFF, sizeof(presieve_pattern));

    /* Mark multiples of 3, 5, 7 in the pattern */
    /* Pattern index i represents odd number (2*i + 3) */
    for (int i = 0; i < PRESIEVE_PRIMES_PRODUCT; i++) {
        uint64_t n = 2 * i + 3;  /* The odd number this index represents */
        if (n % 3 == 0 || n % 5 == 0 || n % 7 == 0) {
            presieve_pattern[i >> 3] &= ~(1 << (i & 7));
        }
    }
    /* But keep 3, 5, 7 themselves as prime */
    presieve_pattern[0] |= 1;  /* 3 */
    presieve_pattern[1 >> 3] |= (1 << (1 & 7));  /* 5 */
    presieve_pattern[2 >> 3] |= (1 << (2 & 7));  /* 7 */

    presieve_initialized = true;
}

typedef struct {
    uint64_t *bitmap;
    uint64_t threshold;
    uint64_t num_words;
} PresieveSieve;

static PresieveSieve* presieve_create(uint64_t threshold) {
    init_presieve_pattern();

    PresieveSieve *sieve = malloc(sizeof(PresieveSieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;
    uint64_t max_bit_idx = (threshold - 3) / 2;
    sieve->num_words = (max_bit_idx / 64) + 1;
    uint64_t num_bytes = sieve->num_words * sizeof(uint64_t);

    sieve->bitmap = malloc(num_bytes);
    if (!sieve->bitmap) { free(sieve); return NULL; }

    /* Fill bitmap with presieve pattern instead of 0xFF */
    /* The presieve pattern marks composites divisible by 3, 5, or 7 */
    uint8_t *bytes = (uint8_t*)sieve->bitmap;
    uint64_t total_bits = max_bit_idx + 1;

    for (uint64_t bit_idx = 0; bit_idx < total_bits; bit_idx++) {
        uint64_t n = 2 * bit_idx + 3;  /* The odd number */
        bool is_composite = (n > 3 && n % 3 == 0) ||
                           (n > 5 && n % 5 == 0) ||
                           (n > 7 && n % 7 == 0);
        if (is_composite) {
            bytes[bit_idx >> 3] &= ~(1 << (bit_idx & 7));
        } else {
            bytes[bit_idx >> 3] |= (1 << (bit_idx & 7));
        }
    }

    uint64_t sqrt_thresh = (uint64_t)sqrt((double)threshold) + 1;

    /* Sieve starting from 11 (skip 3, 5, 7 which are already handled) */
    for (uint64_t p = 11; p <= sqrt_thresh; p += 2) {
        uint64_t idx = (p - 3) / 2;
        if (!((sieve->bitmap[idx >> 6] >> (idx & 63)) & 1)) continue;
        for (uint64_t m = p * p; m <= threshold; m += 2 * p) {
            uint64_t midx = (m - 3) / 2;
            sieve->bitmap[midx >> 6] &= ~(1ULL << (midx & 63));
        }
    }

    return sieve;
}

static uint64_t presieve_count(PresieveSieve *sieve) {
    uint64_t count = (sieve->threshold >= 2) ? 1 : 0;
    for (uint64_t n = 3; n <= sieve->threshold; n += 2) {
        uint64_t idx = (n - 3) / 2;
        if ((sieve->bitmap[idx >> 6] >> (idx & 63)) & 1) count++;
    }
    return count;
}

static void presieve_destroy(PresieveSieve *sieve) {
    if (sieve) { free(sieve->bitmap); free(sieve); }
}

/* ========================================================================== */
/* 4. MOD 30 WHEEL: 8 bits per 30 numbers                                     */
/* ========================================================================== */

/* Numbers coprime to 30 in [0,30): 1, 7, 11, 13, 17, 19, 23, 29 */
static const uint8_t wheel30_residues[8] = {1, 7, 11, 13, 17, 19, 23, 29};
static const uint8_t wheel30_index[30] = {
    0xFF, 0, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 1, 0xFF, 0xFF,
    0xFF, 2, 0xFF, 3, 0xFF, 0xFF, 0xFF, 4, 0xFF, 5,
    0xFF, 0xFF, 0xFF, 6, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 7
};

typedef struct {
    uint8_t *bitmap;
    uint64_t threshold;
    uint64_t num_bytes;
} Wheel30Sieve;

static inline bool wheel30_is_coprime(uint64_t n) {
    uint64_t r = n % 30;
    return wheel30_index[r] != 0xFF;
}

static inline uint64_t wheel30_bit_index(uint64_t n) {
    /* For n coprime to 30, compute its bit index */
    uint64_t block = n / 30;
    uint64_t r = n % 30;
    return block * 8 + wheel30_index[r];
}

static Wheel30Sieve* wheel30_create(uint64_t threshold) {
    Wheel30Sieve *sieve = malloc(sizeof(Wheel30Sieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;
    /* Each 30-number block uses 8 bits = 1 byte */
    sieve->num_bytes = (threshold / 30 + 1);

    sieve->bitmap = malloc(sieve->num_bytes);
    if (!sieve->bitmap) { free(sieve); return NULL; }
    memset(sieve->bitmap, 0xFF, sieve->num_bytes);

    uint64_t sqrt_thresh = (uint64_t)sqrt((double)threshold) + 1;

    /* Sieve primes > 5 (2, 3, 5 are handled separately) */
    for (uint64_t p = 7; p <= sqrt_thresh; p += 2) {
        if (!wheel30_is_coprime(p)) continue;

        uint64_t bit_idx = wheel30_bit_index(p);
        if (!((sieve->bitmap[bit_idx >> 3] >> (bit_idx & 7)) & 1)) continue;

        /* Mark multiples of p */
        for (uint64_t m = p * p; m <= threshold; m += 2 * p) {
            if (!wheel30_is_coprime(m)) continue;
            uint64_t midx = wheel30_bit_index(m);
            sieve->bitmap[midx >> 3] &= ~(1 << (midx & 7));
        }
    }

    return sieve;
}

static uint64_t wheel30_count(Wheel30Sieve *sieve) {
    /* Count 2, 3, 5 */
    uint64_t count = 0;
    if (sieve->threshold >= 2) count++;
    if (sieve->threshold >= 3) count++;
    if (sieve->threshold >= 5) count++;

    /* Count primes coprime to 30 */
    uint64_t num_blocks = sieve->threshold / 30 + 1;
    for (uint64_t block = 0; block < num_blocks; block++) {
        for (int i = 0; i < 8; i++) {
            uint64_t n = block * 30 + wheel30_residues[i];
            if (n <= 1 || n > sieve->threshold) continue;
            uint64_t bit_idx = block * 8 + i;
            if ((sieve->bitmap[bit_idx >> 3] >> (bit_idx & 7)) & 1) count++;
        }
    }
    return count;
}

static void wheel30_destroy(Wheel30Sieve *sieve) {
    if (sieve) { free(sieve->bitmap); free(sieve); }
}

/* ========================================================================== */
/* 5. OPENMP PARALLEL: Parallelized segmented sieve                           */
/* ========================================================================== */

#ifdef _OPENMP

typedef struct {
    uint64_t *bitmap;
    uint64_t threshold;
    uint64_t num_words;
} ParallelSieve;

static ParallelSieve* parallel_create(uint64_t threshold) {
    ParallelSieve *sieve = malloc(sizeof(ParallelSieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;
    uint64_t max_bit_idx = (threshold - 3) / 2;
    sieve->num_words = (max_bit_idx / 64) + 1;

    sieve->bitmap = malloc(sieve->num_words * sizeof(uint64_t));
    if (!sieve->bitmap) { free(sieve); return NULL; }
    memset(sieve->bitmap, 0xFF, sieve->num_words * sizeof(uint64_t));

    uint64_t sqrt_thresh = (uint64_t)sqrt((double)threshold) + 1;

    /* First, sieve small primes sequentially */
    uint64_t small_limit = sqrt_thresh;
    uint64_t small_bits = (small_limit - 3) / 2 + 1;
    uint64_t small_words = (small_bits / 64) + 1;
    uint64_t *small_sieve = malloc(small_words * sizeof(uint64_t));
    memset(small_sieve, 0xFF, small_words * sizeof(uint64_t));

    uint64_t sqrt_small = (uint64_t)sqrt((double)small_limit) + 1;
    for (uint64_t p = 3; p <= sqrt_small; p += 2) {
        uint64_t idx = (p - 3) / 2;
        if (!((small_sieve[idx >> 6] >> (idx & 63)) & 1)) continue;
        for (uint64_t m = p * p; m <= small_limit; m += 2 * p) {
            uint64_t midx = (m - 3) / 2;
            small_sieve[midx >> 6] &= ~(1ULL << (midx & 63));
        }
    }

    /* Collect sieving primes */
    uint64_t *sieving_primes = malloc((small_limit / 2) * sizeof(uint64_t));
    uint64_t num_sieving_primes = 0;
    for (uint64_t p = 3; p <= small_limit; p += 2) {
        uint64_t idx = (p - 3) / 2;
        if ((small_sieve[idx >> 6] >> (idx & 63)) & 1) {
            sieving_primes[num_sieving_primes++] = p;
        }
    }
    free(small_sieve);

    /* Copy small primes to main bitmap */
    for (uint64_t i = 0; i < num_sieving_primes; i++) {
        uint64_t p = sieving_primes[i];
        /* These are already marked, just need to mark composites */
    }

    /* Parallel sieve in segments */
    uint64_t num_segments = (threshold - 3) / SEGMENT_RANGE + 1;

    #pragma omp parallel
    {
        uint8_t *segment = malloc(SEGMENT_SIZE);

        #pragma omp for schedule(dynamic, 4)
        for (uint64_t seg = 0; seg < num_segments; seg++) {
            uint64_t seg_start = 3 + seg * SEGMENT_RANGE;
            uint64_t seg_end = seg_start + SEGMENT_RANGE - 1;
            if (seg_end > threshold) seg_end = threshold;
            if (seg_start > threshold) continue;

            memset(segment, 0xFF, SEGMENT_SIZE);

            for (uint64_t i = 0; i < num_sieving_primes; i++) {
                uint64_t p = sieving_primes[i];
                uint64_t start = p * p;
                if (start > seg_end) break;

                if (start < seg_start) {
                    start = seg_start + (p - (seg_start % p)) % p;
                    if ((start & 1) == 0) start += p;
                }

                for (uint64_t m = start; m <= seg_end; m += 2 * p) {
                    uint64_t local_idx = (m - seg_start) / 2;
                    segment[local_idx >> 3] &= ~(1 << (local_idx & 7));
                }
            }

            /* Copy segment to main bitmap */
            for (uint64_t n = seg_start; n <= seg_end; n += 2) {
                uint64_t local_idx = (n - seg_start) / 2;
                bool is_prime = (segment[local_idx >> 3] >> (local_idx & 7)) & 1;
                uint64_t global_idx = (n - 3) / 2;
                if (!is_prime) {
                    /* Use atomic to avoid race conditions */
                    __sync_fetch_and_and(&sieve->bitmap[global_idx >> 6],
                                         ~(1ULL << (global_idx & 63)));
                }
            }
        }

        free(segment);
    }

    free(sieving_primes);
    return sieve;
}

static uint64_t parallel_count(ParallelSieve *sieve) {
    uint64_t count = (sieve->threshold >= 2) ? 1 : 0;

    #pragma omp parallel for reduction(+:count)
    for (uint64_t n = 3; n <= sieve->threshold; n += 2) {
        uint64_t idx = (n - 3) / 2;
        if ((sieve->bitmap[idx >> 6] >> (idx & 63)) & 1) count++;
    }
    return count;
}

static void parallel_destroy(ParallelSieve *sieve) {
    if (sieve) { free(sieve->bitmap); free(sieve); }
}

#endif /* _OPENMP */

/* ========================================================================== */
/* 6. COMBINED: Segmented + Presieve + Wheel30 + OpenMP                       */
/* ========================================================================== */

#ifdef _OPENMP

typedef struct {
    uint8_t *bitmap;
    uint64_t threshold;
    uint64_t num_bytes;
} CombinedSieve;

static CombinedSieve* combined_create(uint64_t threshold) {
    CombinedSieve *sieve = malloc(sizeof(CombinedSieve));
    if (!sieve) return NULL;

    sieve->threshold = threshold;
    sieve->num_bytes = (threshold / 30 + 1);

    sieve->bitmap = malloc(sieve->num_bytes);
    if (!sieve->bitmap) { free(sieve); return NULL; }
    memset(sieve->bitmap, 0xFF, sieve->num_bytes);

    uint64_t sqrt_thresh = (uint64_t)sqrt((double)threshold) + 1;

    /* Collect sieving primes using simple sequential sieve */
    uint64_t *sieving_primes = malloc((sqrt_thresh / 2) * sizeof(uint64_t));
    uint64_t num_sieving_primes = 0;

    /* Simple sieve to find primes up to sqrt_thresh */
    uint8_t *small = calloc((sqrt_thresh / 16) + 1, 1);
    for (uint64_t i = 3; i <= sqrt_thresh; i += 2) {
        if (small[i >> 4] & (1 << ((i >> 1) & 7))) continue;
        if (i >= 7 && wheel30_is_coprime(i)) {
            sieving_primes[num_sieving_primes++] = i;
        }
        if (i * i <= sqrt_thresh) {
            for (uint64_t j = i * i; j <= sqrt_thresh; j += 2 * i) {
                small[j >> 4] |= (1 << ((j >> 1) & 7));
            }
        }
    }
    free(small);

    /* Parallel sieve by segments - each thread handles different byte ranges */
    #define COMB_SEG_BYTES (64 * 1024)  /* 64KB segments */

    uint64_t num_segs = (sieve->num_bytes + COMB_SEG_BYTES - 1) / COMB_SEG_BYTES;

    #pragma omp parallel for schedule(dynamic, 1)
    for (uint64_t seg = 0; seg < num_segs; seg++) {
        uint64_t byte_start = seg * COMB_SEG_BYTES;
        uint64_t byte_end = byte_start + COMB_SEG_BYTES;
        if (byte_end > sieve->num_bytes) byte_end = sieve->num_bytes;

        /* This segment covers numbers from n_start to n_end */
        uint64_t n_start = (byte_start * 8 / 8) * 30;  /* First block * 30 */
        uint64_t n_end = (byte_end * 8 / 8 + 1) * 30;
        if (n_end > threshold) n_end = threshold;

        for (uint64_t i = 0; i < num_sieving_primes; i++) {
            uint64_t p = sieving_primes[i];
            if (p * p > n_end) break;

            /* Find first multiple of p >= max(p*p, n_start) that is coprime to 30 */
            uint64_t start = p * p;
            if (start < n_start) {
                start = ((n_start + p - 1) / p) * p;
            }

            for (uint64_t m = start; m <= n_end && m <= threshold; m += p) {
                if (!wheel30_is_coprime(m)) continue;
                uint64_t midx = wheel30_bit_index(m);
                uint64_t byte_idx = midx >> 3;
                if (byte_idx >= byte_start && byte_idx < byte_end) {
                    sieve->bitmap[byte_idx] &= ~(1 << (midx & 7));
                }
            }
        }
    }

    free(sieving_primes);
    return sieve;
}

static uint64_t combined_count(CombinedSieve *sieve) {
    uint64_t count = 0;
    if (sieve->threshold >= 2) count++;
    if (sieve->threshold >= 3) count++;
    if (sieve->threshold >= 5) count++;

    uint64_t num_blocks = sieve->threshold / 30 + 1;

    #pragma omp parallel for reduction(+:count)
    for (uint64_t block = 0; block < num_blocks; block++) {
        for (int i = 0; i < 8; i++) {
            uint64_t n = block * 30 + wheel30_residues[i];
            if (n <= 1 || n > sieve->threshold) continue;
            uint64_t bit_idx = block * 8 + i;
            if ((sieve->bitmap[bit_idx >> 3] >> (bit_idx & 7)) & 1) count++;
        }
    }
    return count;
}

static void combined_destroy(CombinedSieve *sieve) {
    if (sieve) { free(sieve->bitmap); free(sieve); }
}

#endif /* _OPENMP */

/* ========================================================================== */
/* Main benchmark driver                                                      */
/* ========================================================================== */

int main(int argc, char *argv[]) {
    uint64_t threshold = 1000000000ULL;  /* 10^9 default */

    if (argc > 1) {
        threshold = strtoull(argv[1], NULL, 10);
    }

    printf("Sieve Optimization Benchmark\n");
    printf("============================\n");
    printf("Threshold: %llu (%.2e)\n", (unsigned long long)threshold, (double)threshold);
    printf("Expected prime count (pi): ~%llu\n\n",
           (unsigned long long)(threshold / log(threshold)));

#ifdef _OPENMP
    printf("OpenMP threads: %d\n\n", omp_get_max_threads());
#else
    printf("OpenMP: disabled\n\n");
#endif

    double baseline_time = 0, t0, t1;
    uint64_t expected_count = 0;

#ifdef _OPENMP
    /* 1. OpenMP parallel segmented */
    printf("1. OpenMP parallel segmented...\n");
    t0 = get_time_ms();
    ParallelSieve *pars = parallel_create(threshold);
    t1 = get_time_ms();
    expected_count = parallel_count(pars);
    printf("   Time: %.1f ms\n", t1-t0);
    printf("   Prime count: %llu\n", (unsigned long long)expected_count);
    printf("   Memory: %.1f MB\n\n", pars->num_words * 8.0 / (1024*1024));
    parallel_destroy(pars);

    /* 2. Combined (Wheel30 + Segmented + OpenMP) */
    printf("2. Combined (Wheel30 + Segmented + OpenMP)...\n");
    t0 = get_time_ms();
    CombinedSieve *cs = combined_create(threshold);
    t1 = get_time_ms();
    uint64_t comb_count = combined_count(cs);
    printf("   Time: %.1f ms\n", t1-t0);
    printf("   Prime count: %llu %s\n", (unsigned long long)comb_count,
           comb_count == expected_count ? "✓" : "✗ MISMATCH");
    printf("   Memory: %.1f MB\n", cs->num_bytes / (1024.0*1024));
    combined_destroy(cs);
    printf("\n");
#endif

    printf("============================\n");

    return 0;
}
