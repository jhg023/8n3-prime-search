/*
 * Optimization 2: Incremental residue updates for trial division
 *
 * Instead of computing candidate % prime for each trial prime,
 * maintain residues incrementally as we iterate through a values.
 *
 * Key insight: As a decreases by 2, candidate increases by delta = 2*(a-1)
 * So: new_candidate mod q = (old_candidate mod q + delta mod q) mod q
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

/* ========================================================================== */
/* Original implementation (from prime.h/solve.h)                             */
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
/* Optimized: Incremental residue updates                                     */
/* ========================================================================== */

/* State for incremental trial division */
typedef struct {
    uint8_t pmod[NUM_TRIAL_PRIMES];  /* candidate mod q for each prime q */
    uint8_t dmod[NUM_TRIAL_PRIMES];  /* delta mod q for each prime q */
} TrialDivState;

/* Initialize state for a given candidate and delta */
static inline void trial_div_init(TrialDivState* state, uint64_t candidate, uint64_t delta) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        uint32_t q = TRIAL_PRIMES[i];
        state->pmod[i] = candidate % q;
        state->dmod[i] = delta % q;
    }
}

/* Update state when a decreases by 2 (candidate increases by delta, delta decreases by 4) */
static inline void trial_div_step(TrialDivState* state) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        uint32_t q = TRIAL_PRIMES[i];
        /* pmod[i] = (pmod[i] + dmod[i]) % q */
        uint32_t new_pmod = state->pmod[i] + state->dmod[i];
        if (new_pmod >= q) new_pmod -= q;
        state->pmod[i] = new_pmod;

        /* dmod[i] = (dmod[i] - 4 + q) % q */
        int32_t new_dmod = (int32_t)state->dmod[i] - 4;
        if (new_dmod < 0) new_dmod += q;
        state->dmod[i] = new_dmod;
    }
}

/* Check if any residue is 0 (meaning divisible by that prime) */
/* Returns: 0 = composite, 1 = is small prime, 2 = needs MR */
static inline int trial_div_check_incremental(TrialDivState* state, uint64_t candidate) {
    for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
        if (state->pmod[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_incr(TrialDivState* state, uint64_t candidate) {
    int td = trial_div_check_incremental(state, candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast_local(candidate);
}

static inline uint64_t find_solution_incr(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    /* Initialize residue state */
    TrialDivState state;
    trial_div_init(&state, candidate, delta);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_incr(&state, candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;

        /* Update candidate and delta */
        candidate += delta;
        delta -= 4;
        a -= 2;

        /* Update residue state */
        trial_div_step(&state);
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

    printf("Original (mod per prime): %.2f n/sec (sum=%llu)\n",
           rate, (unsigned long long)sum);
}

void benchmark_incremental(uint64_t n_start) {
    clock_t start = clock();
    uint64_t sum = 0;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        uint64_t a = find_solution_incr(n_start + i, &p);
        sum += a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    double rate = ITERATIONS / elapsed;

    printf("Incremental residues:     %.2f n/sec (sum=%llu)\n",
           rate, (unsigned long long)sum);
}

/* Verify correctness */
void verify(uint64_t n_start, uint64_t count) {
    int mismatches = 0;
    for (uint64_t i = 0; i < count; i++) {
        uint64_t p_orig, p_incr;
        uint64_t a_orig = find_solution_orig(n_start + i, &p_orig);
        uint64_t a_incr = find_solution_incr(n_start + i, &p_incr);

        if (a_orig != a_incr || p_orig != p_incr) {
            mismatches++;
            if (mismatches <= 5) {
                printf("MISMATCH at n=%llu: orig=(%llu,%llu), incr=(%llu,%llu)\n",
                       (unsigned long long)(n_start + i),
                       (unsigned long long)a_orig, (unsigned long long)p_orig,
                       (unsigned long long)a_incr, (unsigned long long)p_incr);
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

    printf("Testing incremental residue trial division\n");
    printf("Iterations: %d\n\n", ITERATIONS);

    printf("=== Verification ===\n");
    verify(n_start, 10000);
    printf("\n");

    printf("=== n = 10^12 ===\n");
    benchmark_original(n_start);
    benchmark_incremental(n_start);

    printf("\n=== n = 10^15 ===\n");
    n_start = 1000000000000000ULL;
    benchmark_original(n_start);
    benchmark_incremental(n_start);

    printf("\n=== n = 2e18 ===\n");
    n_start = 2000000000000000000ULL;
    benchmark_original(n_start);
    benchmark_incremental(n_start);

    return 0;
}
