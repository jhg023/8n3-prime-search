/*
 * Legendre-based skip: Check (N/ℓ) once per N, then use simple modular check
 *
 * Key insight: If (N/ℓ) = -1, then ℓ NEVER divides any p_a, so we can
 * skip checking divisibility by ℓ entirely (for this N).
 *
 * If (N/ℓ) = +1, we need to check if a² ≡ N (mod ℓ).
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
/* Legendre symbol (fast)                                                     */
/* ========================================================================== */

static inline int legendre(uint64_t a, uint32_t p) {
    a %= p;
    if (a == 0) return 0;
    uint64_t result = powmod64(a, (p - 1) / 2, p);
    return (result == 1) ? 1 : -1;
}

/* ========================================================================== */
/* Trial primes with Legendre precomputation                                  */
/* ========================================================================== */

static const uint32_t PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_PRIMES 30

/*
 * Precomputed data for a given N:
 * - For primes where (N/p) = -1: never check divisibility
 * - For primes where (N/p) = +1: check if a² ≡ N (mod p)
 * - N_mod[i] = N mod p[i] for quick checking
 */
typedef struct {
    uint8_t check_mask;     /* Bit i = 1 if we need to check prime i (Legendre = +1) */
    uint8_t N_mod[NUM_PRIMES];  /* N mod p for each prime */
} LegendreData;

static inline void init_legendre_data(uint64_t N, LegendreData* ld) {
    ld->check_mask = 0;

    for (int i = 0; i < 8 && i < NUM_PRIMES; i++) {  /* Only first 8 fit in mask */
        uint32_t p = PRIMES[i];
        ld->N_mod[i] = N % p;

        if (legendre(N, p) == 1) {
            ld->check_mask |= (1 << i);
        }
    }

    /* Store remaining N_mod values */
    for (int i = 8; i < NUM_PRIMES; i++) {
        ld->N_mod[i] = N % PRIMES[i];
    }
}

/*
 * Check if candidate p_a = (N - a²)/2 is divisible by any of our primes.
 * Uses precomputed Legendre data to skip primes where (N/p) = -1.
 *
 * Returns: 0 = composite, 1 = is small prime, 2 = needs MR
 */
static inline int trial_div_legendre(uint64_t candidate, uint64_t a,
                                      const LegendreData* ld) {
    /* First 8 primes with Legendre optimization */
    for (int i = 0; i < 8; i++) {
        uint32_t p = PRIMES[i];

        /* If (N/p) = -1, p never divides any candidate, skip */
        if (!(ld->check_mask & (1 << i))) continue;

        /* Check if candidate is divisible by p */
        if (candidate % p == 0) {
            return (candidate == p) ? 1 : 0;
        }
    }

    /* Remaining primes (no Legendre optimization, just trial division) */
    for (int i = 8; i < NUM_PRIMES; i++) {
        if (candidate % PRIMES[i] == 0) {
            return (candidate == PRIMES[i]) ? 1 : 0;
        }
    }

    return 2;
}

static inline bool is_prime_legendre(uint64_t candidate, uint64_t a,
                                      const LegendreData* ld) {
    int td = trial_div_legendre(candidate, a, ld);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

/* ========================================================================== */
/* Solution finder with Legendre optimization                                  */
/* ========================================================================== */

static uint64_t find_solution_legendre(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute Legendre data */
    LegendreData ld;
    init_legendre_data(N, &ld);

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;

        uint64_t candidate = (N - a_sq) >> 1;

        if (is_prime_legendre(candidate, a, &ld)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }

    return 0;
}

/* ========================================================================== */
/* Baseline                                                                   */
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
/* Benchmark                                                                  */
/* ========================================================================== */

static void benchmark(const char* name, uint64_t (*fn)(uint64_t, uint64_t*),
                      uint64_t start_n, uint64_t count) {
    double t0 = get_time();
    volatile uint64_t sum = 0;

    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t p;
        uint64_t a = fn(n, &p);
        sum += a;
    }

    double elapsed = get_time() - t0;
    printf("  %-30s %10.0f n/sec  (%.3fs)\n", name, (double)count / elapsed, elapsed);
}

int main(int argc, char** argv) {
    printf("Legendre Symbol Optimization\n");
    printf("============================\n");

    uint64_t count = 200000;
    if (argc > 1) count = strtoull(argv[1], NULL, 10);

    /* Verify */
    printf("\nVerifying...\n");
    for (uint64_t n = 1; n <= 10000; n++) {
        uint64_t p1, p2;
        uint64_t a1 = find_solution_baseline(n, &p1);
        uint64_t a2 = find_solution_legendre(n, &p2);

        uint64_t N = 8 * n + 3;
        if (a1*a1 + 2*p1 != N || a2*a2 + 2*p2 != N) {
            printf("ERROR at n=%" PRIu64 "\n", n);
            return 1;
        }
    }
    printf("OK\n");

    uint64_t scales[] = {
        1000000000ULL,
        1000000000000ULL,
        1000000000000000ULL,
        1000000000000000000ULL,
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        printf("\nBenchmark at n = %.2e\n", (double)scales[i]);
        benchmark("Baseline", find_solution_baseline, scales[i], count);
        benchmark("Legendre skip", find_solution_legendre, scales[i], count);
    }

    return 0;
}
