/*
 * Optimization 2 (Alternative): Bitset sieve for first K candidates
 *
 * Instead of testing candidates one by one, precompute a bitmask
 * of survivors (candidates not divisible by any small prime) for
 * the first K candidates, then only run primality tests on survivors.
 *
 * This is "trial division in parallel" using bitwise operations.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/fj64_table.h"

/* Trial primes (same as prime.h) */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_TRIAL_PRIMES 30

#define ITERATIONS 10000000

/* Batch size - tunable */
#define K 64

/* ========================================================================== */
/* Original implementation                                                    */
/* ========================================================================== */

static inline int trial_division_check_orig(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline bool is_prime_fj64_fast_local(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

static inline bool is_candidate_prime_orig(uint64_t candidate) {
    int td = trial_division_check_orig(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast_local(candidate);
}

static inline uint64_t find_solution_orig(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_orig(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;

        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Bitset sieve implementation                                                */
/* ========================================================================== */

/*
 * Build a bitmask of K candidates that survive trial division.
 * Bit i is set if candidate at index i is NOT divisible by any trial prime.
 * (Index 0 = a_max, index 1 = a_max-2, etc.)
 */
static inline uint64_t build_survivor_mask(uint64_t N, uint64_t a_max, int k) {
    uint64_t survivors = (1ULL << k) - 1;  /* Start with all bits set */

    /* For each trial prime, compute which candidates are divisible */
    for (int pi = 0; pi < NUM_TRIAL_PRIMES; pi++) {
        uint32_t q = TRIAL_PRIMES[pi];

        /* Compute N mod q once */
        uint32_t N_mod_q = N % q;

        /* For each candidate index i:
         * candidate_i = (N - a_i^2) / 2 where a_i = a_max - 2*i
         * candidate_i divisible by q iff (N - a_i^2) mod (2*q) == 0
         * i.e., a_i^2 mod q == N mod q
         *
         * But since we're checking divisibility, we just compute directly.
         */
        uint64_t a = a_max;
        uint64_t candidate = (N - a * a) >> 1;
        uint64_t delta = 2 * (a - 1);

        for (int i = 0; i < k && a >= 1; i++) {
            if (candidate >= 2) {
                if (candidate % q == 0) {
                    /* Candidate is divisible by q */
                    if (candidate != q) {
                        /* Composite - clear bit */
                        survivors &= ~(1ULL << i);
                    }
                    /* If candidate == q, it's prime, keep bit set */
                }
            } else {
                /* Candidate < 2, not valid */
                survivors &= ~(1ULL << i);
            }

            candidate += delta;
            delta -= 4;
            a -= 2;
        }
    }

    return survivors;
}

/* Test primality only on survivors */
static inline uint64_t find_solution_sieve(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    /* Count how many candidates we have */
    uint64_t num_candidates = (a_max + 1) / 2;

    /* Process in batches of K */
    uint64_t a_base = a_max;

    while (a_base >= 1) {
        /* Determine batch size */
        int batch_size = (a_base < K * 2) ? (a_base + 1) / 2 : K;
        if (batch_size > 64) batch_size = 64;
        if (batch_size <= 0) break;

        /* Build survivor mask for this batch */
        uint64_t survivors = build_survivor_mask(N, a_base, batch_size);

        /* Test survivors */
        uint64_t a = a_base;
        uint64_t candidate = (N - a * a) >> 1;
        uint64_t delta = 2 * (a - 1);

        for (int i = 0; i < batch_size && a >= 1; i++) {
            if (survivors & (1ULL << i)) {
                /* Survivor - do full primality test */
                if (candidate >= 2) {
                    /* Skip trial division (already done), just check for small primes */
                    bool is_prime = false;
                    if (candidate <= 127) {
                        is_prime = true;
                    } else {
                        is_prime = is_prime_fj64_fast_local(candidate);
                    }

                    if (is_prime) {
                        if (p_out) *p_out = candidate;
                        return a;
                    }
                }
            }

            candidate += delta;
            delta -= 4;
            a -= 2;
        }

        /* Move to next batch */
        a_base -= 2 * batch_size;
    }

    return 0;
}

/* ========================================================================== */
/* Simpler version: process candidates one at a time but batch trial division */
/* ========================================================================== */

/* For a single candidate, check trial division with early exit */
static inline bool passes_trial_division(uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]);
        }
    }
    return true;
}

static inline uint64_t find_solution_simple_batch(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (passes_trial_division(candidate)) {
                /* Passed trial division, now check primality */
                bool is_prime = false;
                if (candidate <= 127) {
                    is_prime = true;
                } else {
                    is_prime = is_prime_fj64_fast_local(candidate);
                }

                if (is_prime) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            }
        }

        if (a < 3) break;

        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Benchmark                                                                  */
/* ========================================================================== */

void benchmark_original(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        uint64_t a = find_solution_orig(n_start + i, &p);
        sum += a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("Original:              %.2f n/sec (sum=%llu)\n",
           rate, (unsigned long long)sum);
}

void benchmark_sieve(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        uint64_t a = find_solution_sieve(n_start + i, &p);
        sum += a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("Bitset sieve (K=%d):   %.2f n/sec (sum=%llu)\n",
           K, rate, (unsigned long long)sum);
}

/* Verify correctness */
void verify(uint64_t n_start, uint64_t count) {
    int mismatches = 0;
    for (uint64_t i = 0; i < count; i++) {
        uint64_t p_orig, p_sieve;
        uint64_t a_orig = find_solution_orig(n_start + i, &p_orig);
        uint64_t a_sieve = find_solution_sieve(n_start + i, &p_sieve);

        if (a_orig != a_sieve) {
            mismatches++;
            if (mismatches <= 5) {
                printf("MISMATCH at n=%llu: orig=(%llu), sieve=(%llu)\n",
                       (unsigned long long)(n_start + i),
                       (unsigned long long)a_orig,
                       (unsigned long long)a_sieve);
            }
        }
    }
    printf("Verification: %d mismatches out of %llu\n", mismatches, (unsigned long long)count);
}

int main(int argc, char **argv) {
    uint64_t n_start = 1000000000000ULL;
    if (argc > 1) {
        n_start = strtoull(argv[1], NULL, 10);
    }

    printf("Testing bitset sieve trial division\n");
    printf("Iterations: %d, Batch size K=%d\n\n", ITERATIONS, K);

    printf("=== Verification ===\n");
    verify(n_start, 10000);
    printf("\n");

    printf("=== n = 10^12 ===\n");
    benchmark_original(n_start);
    benchmark_sieve(n_start);

    printf("\n=== n = 10^15 ===\n");
    n_start = 1000000000000000ULL;
    benchmark_original(n_start);
    benchmark_sieve(n_start);

    printf("\n=== n = 2e18 ===\n");
    n_start = 2000000000000000000ULL;
    benchmark_original(n_start);
    benchmark_sieve(n_start);

    return 0;
}
