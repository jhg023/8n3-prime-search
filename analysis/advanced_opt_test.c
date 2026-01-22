/*
 * Advanced Optimization Tests
 *
 * 1. Prefetch with earlier trigger
 * 2. Skip even a values check (redundant but compiler might miss)
 * 3. Batch N % small_primes computation
 * 4. 2-way interleaved primality testing
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

/* ========================================================================== */
/* Timing                                                                     */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Baseline                                                                   */
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
/* Optimization 1: Aggressive prefetch                                        */
/* ========================================================================== */

static inline bool is_prime_fj64_aggressive_prefetch(uint64_t n) {
    /* Prefetch hash table entry immediately */
    uint32_t hash_idx = fj64_hash(n);
    __builtin_prefetch(&fj64_bases[hash_idx], 0, 3);

    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);

    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;

    return mr_witness_montgomery_cached(n, fj64_bases[hash_idx], n_inv, r_sq);
}

static inline bool is_candidate_prime_prefetch(uint64_t candidate) {
    int td = trial_division_baseline(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_aggressive_prefetch(candidate);
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
/* Optimization 2: Speculative next-candidate prefetch                        */
/* ========================================================================== */

static uint64_t find_solution_speculative(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            /* Speculatively prefetch for next candidate's hash lookup */
            uint64_t next_candidate = candidate + delta - 4;
            if (next_candidate > 127) {
                uint32_t next_hash = fj64_hash(next_candidate);
                __builtin_prefetch(&fj64_bases[next_hash], 0, 1);
            }

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
/* Optimization 3: Unrolled trial division (4x)                               */
/* ========================================================================== */

static inline int trial_division_unroll4(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    if (candidate % 11 == 0) return (candidate == 11) ? 1 : 0;
    if (candidate % 13 == 0) return (candidate == 13) ? 1 : 0;
    if (candidate % 17 == 0) return (candidate == 17) ? 1 : 0;
    if (candidate % 19 == 0) return (candidate == 19) ? 1 : 0;

    /* Unrolled 4x for remaining primes */
    int i = 7;
    while (i + 3 < NUM_TRIAL_PRIMES) {
        if (candidate % TRIAL_PRIMES[i] == 0)
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+1] == 0)
            return (candidate == TRIAL_PRIMES[i+1]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+2] == 0)
            return (candidate == TRIAL_PRIMES[i+2]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+3] == 0)
            return (candidate == TRIAL_PRIMES[i+3]) ? 1 : 0;
        i += 4;
    }
    while (i < NUM_TRIAL_PRIMES) {
        if (candidate % TRIAL_PRIMES[i] == 0)
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        i++;
    }
    return 2;
}

static inline bool is_candidate_prime_unroll4(uint64_t candidate) {
    int td = trial_division_unroll4(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_aggressive_prefetch(candidate);
}

static uint64_t find_solution_unroll4(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_unroll4(candidate)) {
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
/* Optimization 4: Early return for very small candidates                     */
/* ========================================================================== */

/* Small prime lookup table (primes < 256) */
static const uint8_t SMALL_PRIME_BITS[32] = {
    /* Bitmask for primes: bit i set if 2i+1 is prime (odd numbers only) */
    /* 1,3,5,7,... -> bits 0,1,2,3,... */
    0x28, 0x12, 0x21, 0x02, 0x28, 0x06, 0x20, 0x12,  /* 1-63 */
    0x01, 0x0a, 0x04, 0x20, 0x10, 0x02, 0x21, 0x08,  /* 65-127 */
    0x00, 0x48, 0x00, 0x02, 0x08, 0x10, 0x00, 0x08,  /* 129-191 */
    0x20, 0x40, 0x20, 0x00, 0x01, 0x08, 0x04, 0x00,  /* 193-255 */
};

static inline bool is_small_prime(uint64_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;
    if (n < 256) {
        int idx = (int)(n >> 1);
        return (SMALL_PRIME_BITS[idx >> 3] >> (idx & 7)) & 1;
    }
    return false;  /* Not handled by this function */
}

static inline bool is_candidate_prime_early_small(uint64_t candidate) {
    /* Very fast path for small candidates */
    if (candidate < 256) {
        return is_small_prime(candidate);
    }

    int td = trial_division_baseline(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    return is_prime_fj64_aggressive_prefetch(candidate);
}

static uint64_t find_solution_early_small(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_early_small(candidate)) {
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
/* Optimization 5: Combined best ideas                                        */
/* ========================================================================== */

static inline bool is_candidate_prime_best(uint64_t candidate) {
    /* Very fast path for small candidates using bitmask */
    if (candidate < 256) {
        return is_small_prime(candidate);
    }

    /* Trial division with first 7 inlined */
    if (candidate % 3 == 0) return false;
    if (candidate % 5 == 0) return false;
    if (candidate % 7 == 0) return false;
    if (candidate % 11 == 0) return false;
    if (candidate % 13 == 0) return false;
    if (candidate % 17 == 0) return false;
    if (candidate % 19 == 0) return false;

    for (int i = 7; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return false;
        }
    }

    /* Miller-Rabin with prefetch */
    return is_prime_fj64_aggressive_prefetch(candidate);
}

static uint64_t find_solution_best(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_best(candidate)) {
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
    {"baseline",          find_solution_baseline},
    {"prefetch",          find_solution_prefetch},
    {"speculative",       find_solution_speculative},
    {"unroll4+prefetch",  find_solution_unroll4},
    {"early_small",       find_solution_early_small},
    {"best_combined",     find_solution_best},
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

    for (uint64_t i = 0; i < WARMUP_COUNT && i < count; i++) {
        uint64_t p;
        volatile uint64_t a = solver(N, a_max, &p);
        (void)a;

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

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

    printf("Advanced Single-Threaded Optimization Tests\n");
    printf("Iterations per scale: %s\n\n", fmt_num(count));

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

    return 0;
}
