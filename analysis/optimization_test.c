/*
 * Comprehensive Optimization Test Suite
 *
 * Tests various potential optimizations for single-threaded performance:
 * 1. Wheel factorization (skip candidates divisible by 3,5,7)
 * 2. Different trial division prime counts
 * 3. Loop unrolling in trial division
 * 4. Prefetching FJ64 hash table
 * 5. Specialized 32-bit primality test
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

#include "fmt.h"
#include "arith.h"
#include "arith_montgomery.h"
#include "prime.h"
/* fj64_table.h already included by prime.h */

/* ========================================================================== */
/* Timing                                                                     */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Baseline Implementation (Current Code)                                     */
/* ========================================================================== */

static inline int trial_division_baseline(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;

    for (int i = 3; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_baseline(uint64_t candidate) {
    int td = trial_division_baseline(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_baseline(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_baseline(candidate)) {
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
/* Optimization 1: Wheel Factorization (mod 6 wheel)                          */
/* ========================================================================== */

/*
 * For N = 8n + 3, candidate p = (N - a²)/2
 *
 * We want to skip 'a' values where p is guaranteed composite.
 * A mod-6 wheel skips values where p ≡ 0 (mod 2) or p ≡ 0 (mod 3).
 *
 * Since p = (N - a²)/2, p is already guaranteed odd (N is odd, a is odd).
 * For p ≡ 0 (mod 3): N - a² ≡ 0 (mod 6)
 *
 * N = 8n + 3 ≡ 3 (mod 6) for n even, ≡ 5 (mod 6) for n odd
 * Actually: 8n + 3 = 6n + (2n+3). For n=0: 3 mod 6. For n=1: 11 mod 6 = 5, etc.
 *
 * This is complex - let's try a simpler approach: just check a² mod 6
 */

static inline uint64_t find_solution_wheel6(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    /* Precompute N mod 6 */
    int N_mod6 = (int)(N % 6);

    while (1) {
        if (candidate >= 2) {
            /* Skip if p ≡ 0 (mod 3): check if (N - a²) ≡ 0 (mod 6) */
            /* a² mod 6 cycles: 1,1,3,1,1,3 for a=1,3,5,7,9,11... (odd a) */
            int a_mod6 = (int)(a % 6);
            int a_sq_mod6 = (a_mod6 * a_mod6) % 6;
            int diff_mod6 = (N_mod6 - a_sq_mod6 + 6) % 6;

            /* If diff_mod6 == 0, then p = (N-a²)/2 is divisible by 3 */
            if (diff_mod6 != 0 || candidate == 3) {
                if (is_candidate_prime_baseline(candidate)) {
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
/* Optimization 2: Fewer Trial Division Primes (20 instead of 30)             */
/* ========================================================================== */

static const uint32_t TRIAL_PRIMES_20[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73
};

static inline int trial_division_20(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;

    for (int i = 3; i < 20; i++) {
        if (candidate % TRIAL_PRIMES_20[i] == 0) {
            return (candidate == TRIAL_PRIMES_20[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_20(uint64_t candidate) {
    int td = trial_division_20(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 73) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_td20(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_20(candidate)) {
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
/* Optimization 3: More Trial Division Primes (50)                            */
/* ========================================================================== */

static const uint32_t TRIAL_PRIMES_50[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233
};

static inline int trial_division_50(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;

    for (int i = 3; i < 50; i++) {
        if (candidate % TRIAL_PRIMES_50[i] == 0) {
            return (candidate == TRIAL_PRIMES_50[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_50(uint64_t candidate) {
    int td = trial_division_50(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 233) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_td50(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_50(candidate)) {
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
/* Optimization 4: Prefetch FJ64 Hash Table                                   */
/* ========================================================================== */

static inline bool is_prime_fj64_prefetch(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);

    /* Prefetch the hash table entry while computing base-2 test */
    uint32_t hash_idx = fj64_hash(n);
    __builtin_prefetch(&fj64_bases[hash_idx], 0, 3);

    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;

    return mr_witness_montgomery_cached(n, fj64_bases[hash_idx], n_inv, r_sq);
}

static inline bool is_candidate_prime_prefetch(uint64_t candidate) {
    int td = trial_division_baseline(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_prefetch(candidate);
}

static uint64_t find_solution_prefetch(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_prefetch(candidate)) {
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
/* Optimization 5: Inline First 5 Trial Primes                                */
/* ========================================================================== */

static inline int trial_division_inline5(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    if (candidate % 11 == 0) return (candidate == 11) ? 1 : 0;
    if (candidate % 13 == 0) return (candidate == 13) ? 1 : 0;

    for (int i = 5; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_inline5(uint64_t candidate) {
    int td = trial_division_inline5(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_inline5(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_inline5(candidate)) {
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
/* Optimization 6: Combined - Inline5 + Prefetch                              */
/* ========================================================================== */

static inline bool is_candidate_prime_combined(uint64_t candidate) {
    int td = trial_division_inline5(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_prefetch(candidate);
}

static uint64_t find_solution_combined(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_combined(candidate)) {
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
/* Optimization 7: Early termination on small primes in trial division        */
/* ========================================================================== */

/*
 * Since we iterate a from large to small, candidates start small and grow.
 * Early candidates are often small primes themselves. We can check this first.
 */
static inline bool is_candidate_prime_smallfirst(uint64_t candidate) {
    /* Quick check for very small primes */
    if (candidate <= 127) {
        if (candidate < 2) return false;
        if (candidate == 2) return true;
        if ((candidate & 1) == 0) return false;
        for (int i = 0; i < NUM_TRIAL_PRIMES && TRIAL_PRIMES[i] <= candidate; i++) {
            if (candidate == TRIAL_PRIMES[i]) return true;
            if (candidate % TRIAL_PRIMES[i] == 0) return false;
        }
        return true;  /* candidate <= 127 and passed trial division */
    }

    /* Normal path for larger candidates */
    int td = trial_division_baseline(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_smallfirst(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_smallfirst(candidate)) {
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
/* Benchmark Infrastructure                                                   */
/* ========================================================================== */

#define WARMUP_COUNT    100000
#define DEFAULT_COUNT   1000000

typedef uint64_t (*SolverFunc)(uint64_t N, uint64_t a_max, uint64_t* p_out);

typedef struct {
    const char* name;
    SolverFunc func;
} Solver;

static const Solver SOLVERS[] = {
    {"baseline (current)",   find_solution_baseline},
    {"wheel6",              find_solution_wheel6},
    {"TD=20 primes",        find_solution_td20},
    {"TD=50 primes",        find_solution_td50},
    {"prefetch FJ64",       find_solution_prefetch},
    {"inline 5 primes",     find_solution_inline5},
    {"inline5+prefetch",    find_solution_combined},
    {"small-first check",   find_solution_smallfirst},
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
    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    /* Warmup */
    for (uint64_t i = 0; i < WARMUP_COUNT && i < count; i++) {
        uint64_t p;
        volatile uint64_t a = solver(N, a_max, &p);
        (void)a;

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    /* Reset */
    N = 8 * n_start + 3;
    a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

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
    return count / (end - start);
}

int main(int argc, char** argv) {
    setbuf(stdout, NULL);

    uint64_t count = DEFAULT_COUNT;
    if (argc > 1) {
        count = strtoull(argv[1], NULL, 10);
    }

    printf("Single-Threaded Optimization Test Suite\n");
    printf("Iterations per scale: %s\n\n", fmt_num(count));

    printf("Optimizations tested:\n");
    printf("  1. wheel6: Skip a values where p divisible by 3\n");
    printf("  2. TD=20: Fewer trial division primes\n");
    printf("  3. TD=50: More trial division primes\n");
    printf("  4. prefetch: Prefetch FJ64 hash table during base-2 test\n");
    printf("  5. inline5: Inline first 5 trial primes instead of 3\n");
    printf("  6. combined: inline5 + prefetch\n");
    printf("  7. small-first: Optimize for small candidate primes\n");
    printf("\n");

    /* Header */
    printf("%-8s", "Scale");
    for (size_t s = 0; s < NUM_SOLVERS; s++) {
        printf(" %16s", SOLVERS[s].name);
    }
    printf("\n");
    for (size_t i = 0; i < 8 + NUM_SOLVERS * 17; i++) printf("-");
    printf("\n");

    /* Run benchmarks */
    for (size_t i = 0; i < NUM_SCALES; i++) {
        printf("%-8s", SCALES[i].label);

        double baseline_rate = 0;
        for (size_t s = 0; s < NUM_SOLVERS; s++) {
            double rate = run_benchmark(SOLVERS[s].func, SCALES[i].n_start, count);
            if (s == 0) baseline_rate = rate;

            double speedup = rate / baseline_rate;
            if (speedup >= 1.0) {
                printf(" %8s(+%4.1f%%)", fmt_num((uint64_t)rate), (speedup - 1) * 100);
            } else {
                printf(" %8s(%5.1f%%)", fmt_num((uint64_t)rate), (speedup - 1) * 100);
            }
        }
        printf("\n");
    }

    printf("\nLegend: rate (percent change from baseline)\n");
    printf("Positive = faster, Negative = slower\n");

    return 0;
}
