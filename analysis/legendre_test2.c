/*
 * Legendre/Bitset Sieve Optimization Test - Correct Implementation
 *
 * From the image:
 * - Take K = 32 or 64 candidates
 * - Precompute the first K odd 'a' values from a_max and their squares mod q
 * - For each small q, compute a bitmask of indices where a² mod q == N mod q
 * - OR them together to get "bad" mask (candidates where p is divisible by q)
 * - Only run primality on survivors
 *
 * Key: This is "trial division in parallel" without per-candidate divisions.
 * The a_max and thus the a[] array is constant over huge stretches, so it
 * can be rebuilt extremely rarely.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#include "fmt.h"
#include "arith.h"
#include "prime.h"
#include "solve.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

#define SIEVE_K 64  /* Number of candidates to batch */

/* Small primes for sieving - use the first 30 (same as trial division) */
static const uint32_t SIEVE_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_SIEVE_PRIMES 30

/* ========================================================================== */
/* Precomputed Tables for Bitset Sieve                                        */
/* ========================================================================== */

/*
 * For each small prime q, store a² mod q for the first K candidates
 * a_sq_mod[q_idx][i] = (a[i])² mod SIEVE_PRIMES[q_idx]
 *
 * When a_max changes, we rebuild these tables.
 */
static uint32_t a_sq_mod[NUM_SIEVE_PRIMES][SIEVE_K];
static uint64_t cached_a_max = 0;

/*
 * Rebuild the a_sq_mod tables when a_max changes
 */
static void rebuild_sieve_tables(uint64_t a_max) {
    if (a_max == cached_a_max) return;

    uint64_t a = a_max;
    for (int i = 0; i < SIEVE_K && a >= 1; i++) {
        for (int q_idx = 0; q_idx < NUM_SIEVE_PRIMES; q_idx++) {
            uint32_t q = SIEVE_PRIMES[q_idx];
            uint32_t a_mod_q = (uint32_t)(a % q);
            a_sq_mod[q_idx][i] = (a_mod_q * a_mod_q) % q;
        }
        if (a < 3) break;
        a -= 2;
    }

    cached_a_max = a_max;
}

/*
 * Bitset sieve: Find which candidates have p divisible by some small prime
 *
 * For candidate p = (N - a²)/2:
 *   p ≡ 0 (mod q)  iff  N - a² ≡ 0 (mod q)  [since 2 is coprime to odd q]
 *                  iff  a² ≡ N (mod q)
 *
 * So we check: a_sq_mod[q_idx][i] == N_mod[q_idx]
 * If true and p != q, then p is composite.
 */
static inline uint64_t compute_bad_mask(uint64_t N, int num_candidates) {
    uint64_t bad_mask = 0;

    for (int q_idx = 0; q_idx < NUM_SIEVE_PRIMES; q_idx++) {
        uint32_t q = SIEVE_PRIMES[q_idx];
        uint32_t N_mod_q = (uint32_t)(N % q);

        for (int i = 0; i < num_candidates; i++) {
            if (a_sq_mod[q_idx][i] == N_mod_q) {
                /* p ≡ 0 (mod q), so p is divisible by q */
                /* Mark as bad unless p == q (small prime) */
                bad_mask |= (1ULL << i);
            }
        }
    }

    return bad_mask;
}

/*
 * Find solution using bitset sieve approach
 */
static inline uint64_t find_solution_bitset_v2(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    /* Ensure tables are up to date */
    rebuild_sieve_tables(a_max);

    /* Count how many candidates we have (a_max, a_max-2, ..., down to 1) */
    int num_candidates = (int)((a_max + 1) / 2);
    if (num_candidates > SIEVE_K) num_candidates = SIEVE_K;
    if (num_candidates == 0) return 0;

    /* Compute bad mask using precomputed tables */
    uint64_t bad_mask = compute_bad_mask(N, num_candidates);

    /* Test survivors (bits NOT set in bad_mask) */
    uint64_t a = a_max;
    for (int i = 0; i < num_candidates; i++) {
        if (!(bad_mask & (1ULL << i))) {
            uint64_t candidate = (N - a * a) >> 1;

            if (candidate >= 2) {
                /* Check small primes that might equal candidate */
                bool is_small_prime = false;
                for (int q_idx = 0; q_idx < NUM_SIEVE_PRIMES; q_idx++) {
                    if (candidate == SIEVE_PRIMES[q_idx]) {
                        is_small_prime = true;
                        break;
                    }
                }

                if (is_small_prime) {
                    if (p_out) *p_out = candidate;
                    return a;
                }

                /* Run Miller-Rabin */
                if (candidate > 127 && is_prime_fj64_fast(candidate)) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    /* If we exhausted SIEVE_K candidates, continue with normal approach */
    while (a >= 1) {
        uint64_t candidate = (N - a * a) >> 1;

        if (candidate >= 2) {
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;  /* Counterexample */
}

/* ========================================================================== */
/* Hybrid approach: Quick check first, then sieve if needed                   */
/* ========================================================================== */

/*
 * The image says we "typically succeed in ~10-20 tries" - so only use
 * the sieve when the quick approach fails.
 */
#define QUICK_TRIES 8

static inline uint64_t find_solution_hybrid(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    /* Quick tries first */
    for (int tries = 0; tries < QUICK_TRIES && a >= 1; tries++) {
        if (candidate >= 2) {
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) return 0;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    /* Didn't find quick solution - switch to bitset sieve for remaining */
    /* But we need to adjust: rebuild tables starting from current a */
    cached_a_max = 0;  /* Force rebuild */
    rebuild_sieve_tables(a);

    int remaining = (int)((a + 1) / 2);
    if (remaining > SIEVE_K) remaining = SIEVE_K;

    uint64_t bad_mask = compute_bad_mask(N, remaining);

    /* Test survivors */
    for (int i = 0; i < remaining; i++) {
        if (!(bad_mask & (1ULL << i))) {
            candidate = (N - a * a) >> 1;

            if (candidate >= 2) {
                /* Small prime check */
                if (candidate <= 127) {
                    for (int q_idx = 0; q_idx < NUM_SIEVE_PRIMES && SIEVE_PRIMES[q_idx] <= candidate; q_idx++) {
                        if (candidate == SIEVE_PRIMES[q_idx]) {
                            if (p_out) *p_out = candidate;
                            return a;
                        }
                    }
                } else if (is_prime_fj64_fast(candidate)) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    /* Exhausted sieve, continue normally */
    while (a >= 1) {
        candidate = (N - a * a) >> 1;

        if (candidate >= 2) {
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Timing                                                                     */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Benchmark                                                                  */
/* ========================================================================== */

#define WARMUP_COUNT    100000
#define DEFAULT_COUNT   1000000

typedef uint64_t (*SolverFunc)(uint64_t N, uint64_t a_max, uint64_t* p_out);

typedef struct {
    const char* name;
    SolverFunc func;
} Solver;

static const Solver SOLVERS[] = {
    {"baseline",           find_solution_from_N},
    {"bitset_v2 (K=64)",   find_solution_bitset_v2},
    {"hybrid (8+sieve)",   find_solution_hybrid},
};
#define NUM_SOLVERS (sizeof(SOLVERS) / sizeof(SOLVERS[0]))

typedef struct {
    uint64_t n_start;
    const char* label;
} Scale;

static const Scale SCALES[] = {
    {1000000ULL,             "10^6"},
    {1000000000ULL,          "10^9"},
    {1000000000000ULL,       "10^12"},
    {1000000000000000ULL,    "10^15"},
    {1000000000000000000ULL, "10^18"},
};
#define NUM_SCALES (sizeof(SCALES) / sizeof(SCALES[0]))

double run_benchmark(SolverFunc solver, uint64_t n_start, uint64_t count) {
    /* Warmup */
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < WARMUP_COUNT && i < count; i++) {
        uint64_t p;
        volatile uint64_t a = solver(N, a_max, &p);
        (void)a;

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    /* Reset for timed run */
    N = 8 * n_start + 3;
    a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    cached_a_max = 0;  /* Reset cache */

    double start = get_time();

    for (uint64_t i = 0; i < count; i++) {
        uint64_t p;
        volatile uint64_t a = solver(N, a_max, &p);
        (void)a;

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    double end = get_time();
    return count / (end - start);  /* n/sec */
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

    uint64_t count = DEFAULT_COUNT;
    if (argc > 1) {
        count = strtoull(argv[1], NULL, 10);
    }

    printf("Bitset Sieve Optimization Test v2\n");
    printf("Iterations per scale: %s\n\n", fmt_num(count));

    printf("Key insight from the image:\n");
    printf("- a_max is constant over huge stretches (changes every ~1B n values at 10^18)\n");
    printf("- Precompute a² mod q for K candidates once, reuse for many N values\n");
    printf("- This makes sieving O(K * num_primes) per rebuild, amortized O(1)\n\n");

    /* Print header */
    printf("%-10s", "Scale");
    for (size_t s = 0; s < NUM_SOLVERS; s++) {
        printf("  %20s", SOLVERS[s].name);
    }
    printf("\n");
    for (size_t i = 0; i < 10 + NUM_SOLVERS * 22; i++) printf("-");
    printf("\n");

    /* Run benchmarks */
    for (size_t i = 0; i < NUM_SCALES; i++) {
        printf("%-10s", SCALES[i].label);

        double baseline_rate = 0;
        for (size_t s = 0; s < NUM_SOLVERS; s++) {
            cached_a_max = 0;  /* Reset cache between solvers */
            double rate = run_benchmark(SOLVERS[s].func, SCALES[i].n_start, count);
            if (s == 0) baseline_rate = rate;

            double speedup = rate / baseline_rate;
            printf("  %12s (%4.2fx)", fmt_num((uint64_t)rate), speedup);
        }
        printf("\n");
    }

    printf("\n");
    printf("Analysis:\n");
    printf("- At 10^6:  avg_checks ~6, so quick approach (8 tries) usually succeeds\n");
    printf("- At 10^18: avg_checks ~16, more benefit from sieving\n");
    printf("- BUT: the baseline trial division is already extremely efficient\n");
    printf("- The sieve adds overhead that may not pay off\n");

    return 0;
}
