/*
 * Optimization: Sieve on a-values
 *
 * Instead of checking each a one-by-one, create a sieve:
 * 1. For each small prime ℓ, compute Legendre symbol (N/ℓ)
 * 2. If +1: find roots r₁, r₂ where r² ≡ N (mod ℓ)
 * 3. Mark all a ≡ r₁ or r₂ (mod ℓ) as eliminated
 * 4. Only test primality for unmarked a's
 *
 * Use segmented sieve to handle large √N without huge memory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <inttypes.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/prime.h"
#include "../include/solve.h"

/* ========================================================================== */
/* Timing                                                                     */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Legendre symbol and Tonelli-Shanks                                         */
/* ========================================================================== */

/* Legendre symbol (a/p): returns 1 if QR, -1 if NQR, 0 if a ≡ 0 */
static inline int legendre(uint64_t a, uint64_t p) {
    a %= p;
    if (a == 0) return 0;
    uint64_t result = powmod64(a, (p - 1) / 2, p);
    return (result == 1) ? 1 : -1;
}

/* Tonelli-Shanks: find r such that r² ≡ n (mod p) */
static uint64_t tonelli_shanks(uint64_t n, uint64_t p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    n %= p;
    if (n == 0) return 0;

    /* Check if n is a quadratic residue */
    if (powmod64(n, (p - 1) / 2, p) != 1) {
        return 0;  /* No square root exists */
    }

    /* Special case: p ≡ 3 (mod 4) */
    if ((p & 3) == 3) {
        return powmod64(n, (p + 1) / 4, p);
    }

    /* General Tonelli-Shanks */
    uint64_t Q = p - 1;
    int S = 0;
    while ((Q & 1) == 0) {
        Q >>= 1;
        S++;
    }

    /* Find quadratic non-residue */
    uint64_t z = 2;
    while (powmod64(z, (p - 1) / 2, p) != p - 1) {
        z++;
    }

    int M = S;
    uint64_t c = powmod64(z, Q, p);
    uint64_t t = powmod64(n, Q, p);
    uint64_t R = powmod64(n, (Q + 1) / 2, p);

    while (1) {
        if (t == 1) return R;

        int i = 1;
        uint64_t temp = mulmod64(t, t, p);
        while (temp != 1) {
            temp = mulmod64(temp, temp, p);
            i++;
        }

        uint64_t b = c;
        for (int j = 0; j < M - i - 1; j++) {
            b = mulmod64(b, b, p);
        }
        M = i;
        c = mulmod64(b, b, p);
        t = mulmod64(t, c, p);
        R = mulmod64(R, b, p);
    }
}

/* ========================================================================== */
/* Sieve primes and roots                                                     */
/* ========================================================================== */

/* Primes for sieving - chosen so their product is manageable */
#define MAX_SIEVE_PRIMES 50
static const uint32_t SIEVE_PRIMES[MAX_SIEVE_PRIMES] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233
};

/* Precomputed sieve data for a given N */
typedef struct {
    int num_primes;           /* Number of primes with (N/ℓ) = +1 */
    uint32_t primes[MAX_SIEVE_PRIMES];
    uint32_t root1[MAX_SIEVE_PRIMES];  /* First root r₁ */
    uint32_t root2[MAX_SIEVE_PRIMES];  /* Second root r₂ = ℓ - r₁ */
} SieveData;

static void init_sieve_data(uint64_t N, SieveData* sd, int num_sieve_primes) {
    sd->num_primes = 0;

    for (int i = 0; i < num_sieve_primes && i < MAX_SIEVE_PRIMES; i++) {
        uint32_t ell = SIEVE_PRIMES[i];
        int leg = legendre(N, ell);

        if (leg == 1) {
            /* N is a QR mod ℓ, find the roots */
            uint64_t r = tonelli_shanks(N % ell, ell);
            sd->primes[sd->num_primes] = ell;
            sd->root1[sd->num_primes] = (uint32_t)r;
            sd->root2[sd->num_primes] = (uint32_t)(ell - r);
            sd->num_primes++;
        }
        /* If leg == -1, ℓ never divides any p_a, skip it */
        /* If leg == 0, N ≡ 0 (mod ℓ), handle specially if needed */
    }
}

/* ========================================================================== */
/* Segmented Sieve Implementation                                             */
/* ========================================================================== */

#define SEGMENT_SIZE 4096  /* Size of each sieve segment (in a-values) */

/*
 * Sieve a segment of a-values from a_start down to a_start - SEGMENT_SIZE + 1
 * Mark sieve[i] = false if a = a_start - i should be skipped
 */
static void sieve_segment(bool* sieve, uint64_t a_start, int segment_size,
                          const SieveData* sd) {
    /* Initialize all as valid */
    memset(sieve, true, segment_size);

    /* For each sieving prime */
    for (int i = 0; i < sd->num_primes; i++) {
        uint32_t ell = sd->primes[i];
        uint32_t r1 = sd->root1[i];
        uint32_t r2 = sd->root2[i];

        /* Find first a in segment where a ≡ r1 (mod ℓ) */
        /* a_start - offset ≡ r1 (mod ℓ) => offset ≡ a_start - r1 (mod ℓ) */
        uint64_t a_mod = a_start % ell;

        /* For r1: offset where a_start - offset ≡ r1 (mod ℓ) */
        int64_t offset1 = (int64_t)((a_mod - r1 + ell) % ell);
        for (int64_t j = offset1; j < segment_size; j += ell) {
            sieve[j] = false;
        }

        /* For r2 (if different from r1) */
        if (r1 != r2) {
            int64_t offset2 = (int64_t)((a_mod - r2 + ell) % ell);
            for (int64_t j = offset2; j < segment_size; j += ell) {
                sieve[j] = false;
            }
        }
    }
}

/* ========================================================================== */
/* Solution finder with sieve                                                 */
/* ========================================================================== */

static uint64_t find_solution_sieve(uint64_t n, uint64_t* p_out, int num_sieve_primes) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Ensure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    /* Initialize sieve data for this N */
    SieveData sd;
    init_sieve_data(N, &sd, num_sieve_primes);

    /* Sieve segment (only odd a's, so effective segment is 2x) */
    bool sieve[SEGMENT_SIZE];
    uint64_t segment_start = a_max;

    while (segment_start >= 1) {
        /* Calculate segment bounds */
        int segment_size = (segment_start >= SEGMENT_SIZE) ? SEGMENT_SIZE : (int)(segment_start + 1);

        /* Sieve this segment */
        sieve_segment(sieve, segment_start, segment_size, &sd);

        /* Iterate through unmarked a's in this segment (reverse order) */
        for (int i = 0; i < segment_size; i++) {
            uint64_t a = segment_start - i;

            /* Skip even a's */
            if ((a & 1) == 0) continue;

            /* Skip sieved a's (but handle edge case where p_a = ℓ) */
            if (!sieve[i]) {
                /* Check if candidate equals one of the sieve primes */
                uint64_t a_sq = a * a;
                if (a_sq <= N - 4) {
                    uint64_t candidate = (N - a_sq) >> 1;
                    /* Only check if candidate is a small prime */
                    if (candidate <= SIEVE_PRIMES[num_sieve_primes - 1]) {
                        for (int j = 0; j < num_sieve_primes; j++) {
                            if (candidate == SIEVE_PRIMES[j]) {
                                if (p_out) *p_out = candidate;
                                return a;
                            }
                        }
                    }
                }
                continue;
            }

            /* Test this a */
            uint64_t a_sq = a * a;
            if (a_sq <= N - 4) {
                uint64_t candidate = (N - a_sq) >> 1;

                if (is_candidate_prime(candidate)) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            }
        }

        if (segment_start < SEGMENT_SIZE) break;
        segment_start -= SEGMENT_SIZE;
    }

    return 0;  /* Counterexample */
}

/* ========================================================================== */
/* Baseline for comparison                                                    */
/* ========================================================================== */

static uint64_t find_solution_baseline(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;

        uint64_t candidate = (N - a_sq) >> 1;

        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }

    return 0;
}

/* ========================================================================== */
/* Benchmarking                                                               */
/* ========================================================================== */

static void benchmark(const char* name, uint64_t (*fn)(uint64_t, uint64_t*),
                      uint64_t start_n, uint64_t count) {
    double t0 = get_time();
    volatile uint64_t sum = 0;

    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t p;
        uint64_t a = fn(n, &p);
        sum += a;

        /* Verify */
        if (a > 0) {
            uint64_t N = 8 * n + 3;
            if (a * a + 2 * p != N) {
                printf("ERROR in %s at n=%" PRIu64 "\n", name, n);
                exit(1);
            }
        }
    }

    double elapsed = get_time() - t0;
    printf("  %-30s %10.0f n/sec  (%.3fs)\n", name, (double)count / elapsed, elapsed);
}

/* Wrapper functions for different sieve prime counts */
static uint64_t find_solution_sieve_10(uint64_t n, uint64_t* p_out) {
    return find_solution_sieve(n, p_out, 10);
}

static uint64_t find_solution_sieve_20(uint64_t n, uint64_t* p_out) {
    return find_solution_sieve(n, p_out, 20);
}

static uint64_t find_solution_sieve_30(uint64_t n, uint64_t* p_out) {
    return find_solution_sieve(n, p_out, 30);
}

static uint64_t find_solution_sieve_50(uint64_t n, uint64_t* p_out) {
    return find_solution_sieve(n, p_out, 50);
}

int main(int argc, char** argv) {
    printf("Sieve on a-values Optimization\n");
    printf("==============================\n");

    uint64_t count = 100000;
    if (argc > 1) count = strtoull(argv[1], NULL, 10);

    /* Verify correctness */
    printf("\nVerifying correctness (n = 1 to 10000)...\n");
    for (uint64_t n = 1; n <= 10000; n++) {
        uint64_t p1, p2;
        uint64_t a1 = find_solution_baseline(n, &p1);
        uint64_t a2 = find_solution_sieve(n, &p2, 30);

        uint64_t N = 8 * n + 3;
        if (a1 == 0 || a2 == 0) {
            printf("ERROR: No solution for n=%" PRIu64 "\n", n);
            return 1;
        }
        if (a1*a1 + 2*p1 != N || a2*a2 + 2*p2 != N) {
            printf("ERROR: Invalid solution for n=%" PRIu64 "\n", n);
            return 1;
        }
    }
    printf("Verification passed.\n");

    /* Benchmark at multiple scales */
    uint64_t scales[] = {
        1000000000ULL,
        1000000000000ULL,
        1000000000000000ULL,
        1000000000000000000ULL,
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        printf("\nBenchmark at n = %.2e, count = %" PRIu64 "\n",
               (double)scales[i], count);
        printf("%-32s %15s\n", "Implementation", "Throughput");
        printf("--------------------------------------------------------\n");

        benchmark("Baseline", find_solution_baseline, scales[i], count);
        benchmark("Sieve (10 primes)", find_solution_sieve_10, scales[i], count);
        benchmark("Sieve (20 primes)", find_solution_sieve_20, scales[i], count);
        benchmark("Sieve (30 primes)", find_solution_sieve_30, scales[i], count);
        benchmark("Sieve (50 primes)", find_solution_sieve_50, scales[i], count);
    }

    return 0;
}
