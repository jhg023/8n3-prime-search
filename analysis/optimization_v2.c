/*
 * Optimization v2: Direct divisibility check instead of Tonelli-Shanks
 *
 * Key insight: Instead of computing sqrt(N mod q) via Tonelli-Shanks,
 * we can directly check if (a² - N) ≡ 0 (mod q), which means candidate
 * p = (N - a²)/2 is divisible by q.
 *
 * This is O(1) per prime instead of O(log q) for Tonelli-Shanks.
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
/* Timing utilities                                                           */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Strategy 0: Baseline (current implementation)                              */
/* ========================================================================== */

static uint64_t find_solution_baseline(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;

        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;

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
/* Strategy 1: Inline divisibility sieve                                      */
/* ========================================================================== */

/*
 * Instead of trial division on candidate, check if a produces a candidate
 * divisible by small primes. This is done by checking (a² - N) % (2q) == 0.
 *
 * If (N - a²) is divisible by 2q, then (N - a²)/2 is divisible by q.
 */

/* Small primes for inline sieve (avoid those already in trial division) */
#define INLINE_SIEVE_PRIMES_COUNT 8
static const uint32_t INLINE_SIEVE_PRIMES[INLINE_SIEVE_PRIMES_COUNT] = {
    3, 5, 7, 11, 13, 17, 19, 23
};

static inline bool quick_composite_check(uint64_t N, uint64_t a_sq) {
    /* Check if (N - a_sq)/2 is divisible by small primes */
    uint64_t diff = N - a_sq;  /* Already know a_sq <= N - 4 */

    for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT; i++) {
        uint32_t q = INLINE_SIEVE_PRIMES[i];
        uint32_t q2 = q * 2;
        /* diff is even (N is odd, a_sq is odd), so diff/2 is integer */
        /* diff/2 % q == 0 iff diff % (2q) == 0 */
        if (diff % q2 == 0) {
            /* candidate is divisible by q. Is it equal to q? */
            uint64_t candidate = diff >> 1;
            if (candidate != q) {
                return true;  /* Definitely composite */
            }
        }
    }
    return false;  /* Might be prime */
}

static uint64_t find_solution_inline_sieve(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;

        if (a_sq <= N - 4) {
            /* Quick check before full primality test */
            if (!quick_composite_check(N, a_sq)) {
                uint64_t candidate = (N - a_sq) >> 1;

                if (is_candidate_prime(candidate)) {
                    if (p_out) *p_out = candidate;
                    return a;
                }
            }
        }

        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Strategy 2: Precompute N mod 2q for faster modulo                          */
/* ========================================================================== */

static uint64_t find_solution_precomputed_mod(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute N mod 2q for each sieve prime */
    uint32_t N_mod[INLINE_SIEVE_PRIMES_COUNT];
    for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT; i++) {
        N_mod[i] = N % (INLINE_SIEVE_PRIMES[i] * 2);
    }

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;

        if (a_sq <= N - 4) {
            bool skip = false;
            uint64_t candidate = (N - a_sq) >> 1;

            /* Quick composite check using precomputed moduli */
            for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT && !skip; i++) {
                uint32_t q = INLINE_SIEVE_PRIMES[i];
                uint32_t q2 = q * 2;
                uint32_t a_sq_mod = (a_sq) % q2;
                /* (N - a_sq) % q2 == 0 ? */
                if (a_sq_mod == N_mod[i]) {
                    if (candidate != q) {
                        skip = true;
                    }
                }
            }

            if (!skip && is_candidate_prime(candidate)) {
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
/* Strategy 3: Incremental a² computation                                     */
/* ========================================================================== */

/*
 * Instead of computing a² each iteration, use the identity:
 * (a-2)² = a² - 4a + 4
 *
 * We can incrementally update a² as we iterate.
 */

static uint64_t find_solution_incremental(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    uint64_t a_sq = a * a;

    while (1) {
        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;

            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;

        /* Incremental update: (a-2)² = a² - 4a + 4 = a² - 4(a-1) */
        a_sq = a_sq - 4 * a + 4;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Strategy 4: Combined incremental + inline sieve                            */
/* ========================================================================== */

static uint64_t find_solution_combined_v2(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute N mod 2q */
    uint32_t N_mod[INLINE_SIEVE_PRIMES_COUNT];
    for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT; i++) {
        N_mod[i] = N % (INLINE_SIEVE_PRIMES[i] * 2);
    }

    uint64_t a = a_max;
    uint64_t a_sq = a * a;

    while (1) {
        if (a_sq <= N - 4) {
            bool skip = false;
            uint64_t candidate = (N - a_sq) >> 1;

            /* Quick composite check */
            for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT && !skip; i++) {
                uint32_t q = INLINE_SIEVE_PRIMES[i];
                uint32_t q2 = q * 2;
                uint32_t a_sq_mod = a_sq % q2;
                if (a_sq_mod == N_mod[i] && candidate != q) {
                    skip = true;
                }
            }

            if (!skip && is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        a_sq = a_sq - 4 * a + 4;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Strategy 5: Skip redundant trial division primes                           */
/* ========================================================================== */

/*
 * If we're doing inline sieve for primes 3,5,7,11,13,17,19,23,
 * we can remove those from the trial division list.
 */

static const uint32_t REDUCED_TRIAL_PRIMES[] = {
    29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
    71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
static const int NUM_REDUCED_TRIAL_PRIMES = 22;

static inline int reduced_trial_division(uint64_t candidate) {
    for (int i = 0; i < NUM_REDUCED_TRIAL_PRIMES; i++) {
        if (candidate % REDUCED_TRIAL_PRIMES[i] == 0) {
            return (candidate == REDUCED_TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_reduced(uint64_t candidate) {
    int td = reduced_trial_division(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_reduced_trial(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute N mod 2q */
    uint32_t N_mod[INLINE_SIEVE_PRIMES_COUNT];
    for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT; i++) {
        N_mod[i] = N % (INLINE_SIEVE_PRIMES[i] * 2);
    }

    uint64_t a = a_max;
    uint64_t a_sq = a * a;

    while (1) {
        if (a_sq <= N - 4) {
            bool skip = false;
            uint64_t candidate = (N - a_sq) >> 1;

            /* Inline sieve for small primes */
            for (int i = 0; i < INLINE_SIEVE_PRIMES_COUNT && !skip; i++) {
                uint32_t q = INLINE_SIEVE_PRIMES[i];
                uint32_t q2 = q * 2;
                uint32_t a_sq_mod = a_sq % q2;
                if (a_sq_mod == N_mod[i] && candidate != q) {
                    skip = true;
                }
            }

            if (!skip && is_candidate_prime_reduced(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        a_sq = a_sq - 4 * a + 4;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Benchmarking infrastructure                                                */
/* ========================================================================== */

typedef uint64_t (*solver_fn)(uint64_t n, uint64_t* p_out);

typedef struct {
    const char* name;
    solver_fn solve;
} Strategy;

static Strategy strategies[] = {
    {"Baseline", find_solution_baseline},
    {"Inline sieve (8 primes)", find_solution_inline_sieve},
    {"Precomputed mod", find_solution_precomputed_mod},
    {"Incremental a²", find_solution_incremental},
    {"Combined (incr + sieve)", find_solution_combined_v2},
    {"Reduced trial div", find_solution_reduced_trial},
};

#define NUM_STRATEGIES (sizeof(strategies) / sizeof(strategies[0]))

static void benchmark_strategy(Strategy* s, uint64_t start_n, uint64_t count) {
    double t0 = get_time();
    uint64_t solutions = 0;

    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t p;
        uint64_t a = s->solve(n, &p);
        if (a > 0) {
            solutions++;
            /* Verify */
            uint64_t N = 8 * n + 3;
            if (a * a + 2 * p != N) {
                printf("ERROR: Verification failed for n=%" PRIu64 "\n", n);
                exit(1);
            }
        }
    }

    double elapsed = get_time() - t0;
    double rate = count / elapsed;

    printf("  %-28s %10.0f n/sec  (%.3fs)\n", s->name, rate, elapsed);
}

static void run_benchmarks(uint64_t start_n, uint64_t count) {
    printf("\nBenchmark at n = %.2e, count = %" PRIu64 "\n", (double)start_n, count);
    printf("%-30s %15s\n", "Strategy", "Throughput");
    printf("------------------------------------------------------\n");

    for (size_t i = 0; i < NUM_STRATEGIES; i++) {
        benchmark_strategy(&strategies[i], start_n, count);
    }
}

/* ========================================================================== */
/* Main                                                                       */
/* ========================================================================== */

int main(int argc, char** argv) {
    printf("Optimization v2: Inline Sieve & Incremental Computation\n");
    printf("========================================================\n");

    uint64_t count = 100000;

    if (argc > 1) {
        count = strtoull(argv[1], NULL, 10);
    }

    /* Verify correctness first */
    printf("\nVerifying correctness (n = 1 to 10000)...\n");
    for (uint64_t n = 1; n <= 10000; n++) {
        uint64_t N = 8 * n + 3;

        for (size_t i = 0; i < NUM_STRATEGIES; i++) {
            uint64_t p;
            uint64_t a = strategies[i].solve(n, &p);

            if (a == 0) {
                printf("ERROR: %s found no solution for n=%" PRIu64 "\n",
                       strategies[i].name, n);
                return 1;
            }

            if (a * a + 2 * p != N) {
                printf("ERROR: %s invalid solution for n=%" PRIu64 "\n",
                       strategies[i].name, n);
                return 1;
            }
        }
    }
    printf("All strategies produce valid solutions.\n");

    /* Benchmark at multiple scales */
    uint64_t scales[] = {
        1000000000ULL,           /* 10^9 */
        1000000000000ULL,        /* 10^12 */
        1000000000000000ULL,     /* 10^15 */
        1000000000000000000ULL,  /* 10^18 */
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        run_benchmarks(scales[i], count);
    }

    return 0;
}
