/*
 * Legendre Skip Optimization Test - True Implementation
 *
 * The image suggests: For each small prime q, precompute sqrt[q][r] tables.
 * Then (N/q) = -1 iff sqrt[q][N%q] == NONE.
 *
 * But wait - we're not trying to check if N is a QR. We're trying to find
 * a such that (N - a²)/2 is prime.
 *
 * Let me re-read the image more carefully...
 *
 * The optimization is about checking candidates p = (N - a²)/2.
 * For p to be divisible by q: N - a² ≡ 0 (mod 2q) => a² ≡ N (mod 2q) for odd q
 * Actually for odd q: a² ≡ N (mod q)
 *
 * So for each q, we can precompute which residues r have a² ≡ r (mod q),
 * and then check if N mod q is one of those residues.
 *
 * If a² ≡ N (mod q) for some candidate a, then p = (N-a²)/2 ≡ 0 (mod q),
 * meaning p is divisible by q.
 *
 * The insight: Instead of doing trial division on each p candidate,
 * we can check if a² ≡ N (mod q) to detect divisibility by q.
 *
 * But this is exactly what the current code does implicitly through
 * the trial division `candidate % q == 0`.
 *
 * Let me try a different approach: SKIP candidates entirely based on
 * which 'a' values are "good" vs "bad".
 *
 * An 'a' value is "bad" mod q if a² ≡ N (mod q) AND the resulting p is
 * divisible by q but not equal to q.
 *
 * We can precompute: For each q, which a residues lead to bad candidates?
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
/* Precomputed Square Root Tables                                             */
/* ========================================================================== */

/*
 * Small primes for Legendre checking
 */
static const uint16_t SMALL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_SMALL_PRIMES 30
#define MAX_SMALL_PRIME 127

#define SQRT_NONE 0xFF

/*
 * sqrt_table[q_idx][r] = x where x² ≡ r (mod q), or SQRT_NONE if no such x
 * For q ≤ 127, we need at most 127 bytes per prime = ~4KB total
 */
static uint8_t sqrt_table[NUM_SMALL_PRIMES][MAX_SMALL_PRIME + 1];

static void init_sqrt_tables(void) {
    memset(sqrt_table, SQRT_NONE, sizeof(sqrt_table));

    for (int q_idx = 0; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        uint16_t q = SMALL_PRIMES[q_idx];

        /* 0² ≡ 0 (mod q) */
        sqrt_table[q_idx][0] = 0;

        /* For x in [1, q-1], compute x² mod q and record sqrt */
        for (uint16_t x = 1; x < q; x++) {
            uint16_t sq = (x * x) % q;
            if (sqrt_table[q_idx][sq] == SQRT_NONE) {
                sqrt_table[q_idx][sq] = (uint8_t)x;
            }
        }
    }
}

/*
 * Check if r is a quadratic residue mod SMALL_PRIMES[q_idx]
 */
static inline bool is_qr(int q_idx, uint16_t r) {
    return sqrt_table[q_idx][r] != SQRT_NONE;
}

/*
 * Get sqrt(r) mod SMALL_PRIMES[q_idx], or SQRT_NONE if not a QR
 */
static inline uint8_t get_sqrt(int q_idx, uint16_t r) {
    return sqrt_table[q_idx][r];
}

/* ========================================================================== */
/* Legendre-based candidate checking                                          */
/* ========================================================================== */

/*
 * For candidate p = (N - a²)/2:
 *   p ≡ 0 (mod q) iff a² ≡ N (mod q)
 *
 * So we can check divisibility by q without computing p % q:
 * just check if (a % q)² % q == N % q
 *
 * This trades 1 modulo (p % q) for 2 modulos (a % q, N % q) plus a multiply.
 * That's not a win unless we can precompute N % q.
 *
 * Better: precompute N % q for all q at the start, then for each a:
 *   check if (a % q)² % q == N_mod[q_idx]
 *
 * This is 1 modulo (a % q) plus 1 multiply plus 1 lookup.
 * Original: 1 modulo (candidate % q).
 * Since candidate is larger than a, the Legendre approach might be slightly faster.
 */

/* Precomputed N mod q for each small prime */
static uint16_t N_mod_q[NUM_SMALL_PRIMES];

static void precompute_N_mod(uint64_t N) {
    for (int q_idx = 0; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        N_mod_q[q_idx] = (uint16_t)(N % SMALL_PRIMES[q_idx]);
    }
}

/*
 * Check if p = (N - a²)/2 is divisible by SMALL_PRIMES[q_idx]
 * using precomputed N_mod_q and sqrt tables
 */
static inline bool is_div_by_q(int q_idx, uint64_t a) {
    uint16_t q = SMALL_PRIMES[q_idx];
    uint16_t a_mod_q = (uint16_t)(a % q);
    uint16_t a_sq_mod_q = (a_mod_q * a_mod_q) % q;
    return (a_sq_mod_q == N_mod_q[q_idx]);
}

/*
 * Trial division using Legendre-style checking
 * Returns: 0 = composite, 1 = is small prime q, 2 = needs Miller-Rabin
 */
static inline int trial_division_legendre(uint64_t a, uint64_t candidate) {
    for (int q_idx = 0; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        if (is_div_by_q(q_idx, a)) {
            /* p is divisible by q */
            return (candidate == SMALL_PRIMES[q_idx]) ? 1 : 0;
        }
    }
    return 2;  /* Passed trial division */
}

static inline bool is_candidate_prime_legendre(uint64_t a, uint64_t candidate) {
    int td = trial_division_legendre(a, candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

/*
 * Find solution using Legendre-based trial division
 */
static uint64_t find_solution_legendre(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    /* Precompute N mod q for all small primes */
    precompute_N_mod(N);

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_legendre(a, candidate)) {
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
/* Optimized version: Precompute a % q incrementally                          */
/* ========================================================================== */

/*
 * Since a decreases by 2 each step, we can track a % q incrementally:
 *   new_a_mod = (old_a_mod - 2 + q) % q
 *
 * This avoids the expensive 64-bit modulo operation.
 */

static uint16_t a_mod_q[NUM_SMALL_PRIMES];

static void init_a_mod(uint64_t a_start) {
    for (int q_idx = 0; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        a_mod_q[q_idx] = (uint16_t)(a_start % SMALL_PRIMES[q_idx]);
    }
}

static inline void decrement_a_mod(void) {
    for (int q_idx = 0; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        uint16_t q = SMALL_PRIMES[q_idx];
        /* a -= 2, so a_mod -= 2 (mod q) */
        a_mod_q[q_idx] = (a_mod_q[q_idx] + q - 2) % q;
    }
}

static inline bool is_div_by_q_fast(int q_idx) {
    uint16_t q = SMALL_PRIMES[q_idx];
    uint16_t a_sq_mod_q = (a_mod_q[q_idx] * a_mod_q[q_idx]) % q;
    return (a_sq_mod_q == N_mod_q[q_idx]);
}

static inline int trial_division_legendre_fast(uint64_t candidate) {
    for (int q_idx = 0; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        if (is_div_by_q_fast(q_idx)) {
            return (candidate == SMALL_PRIMES[q_idx]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_legendre_fast(uint64_t candidate) {
    int td = trial_division_legendre_fast(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_legendre_fast(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    precompute_N_mod(N);
    init_a_mod(a_max);

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_legendre_fast(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
        decrement_a_mod();
    }

    return 0;
}

/* ========================================================================== */
/* Even more optimized: Inline the first 3 primes                            */
/* ========================================================================== */

static uint16_t a_mod_3, a_mod_5, a_mod_7;

static void init_a_mod_inline(uint64_t a_start) {
    a_mod_3 = (uint16_t)(a_start % 3);
    a_mod_5 = (uint16_t)(a_start % 5);
    a_mod_7 = (uint16_t)(a_start % 7);
    /* Initialize the rest */
    for (int q_idx = 3; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        a_mod_q[q_idx] = (uint16_t)(a_start % SMALL_PRIMES[q_idx]);
    }
}

static inline void decrement_a_mod_inline(void) {
    a_mod_3 = (a_mod_3 + 1) % 3;  /* -2 ≡ +1 (mod 3) */
    a_mod_5 = (a_mod_5 + 3) % 5;  /* -2 ≡ +3 (mod 5) */
    a_mod_7 = (a_mod_7 + 5) % 7;  /* -2 ≡ +5 (mod 7) */
    for (int q_idx = 3; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        uint16_t q = SMALL_PRIMES[q_idx];
        a_mod_q[q_idx] = (a_mod_q[q_idx] + q - 2) % q;
    }
}

static inline int trial_division_legendre_inline(uint64_t candidate) {
    /* Inline check for q=3: a² ≡ N (mod 3) */
    uint16_t a_sq_3 = (a_mod_3 * a_mod_3) % 3;
    if (a_sq_3 == N_mod_q[0]) {
        return (candidate == 3) ? 1 : 0;
    }

    /* Inline check for q=5 */
    uint16_t a_sq_5 = (a_mod_5 * a_mod_5) % 5;
    if (a_sq_5 == N_mod_q[1]) {
        return (candidate == 5) ? 1 : 0;
    }

    /* Inline check for q=7 */
    uint16_t a_sq_7 = (a_mod_7 * a_mod_7) % 7;
    if (a_sq_7 == N_mod_q[2]) {
        return (candidate == 7) ? 1 : 0;
    }

    /* Loop for remaining primes */
    for (int q_idx = 3; q_idx < NUM_SMALL_PRIMES; q_idx++) {
        uint16_t q = SMALL_PRIMES[q_idx];
        uint16_t a_sq_q = (a_mod_q[q_idx] * a_mod_q[q_idx]) % q;
        if (a_sq_q == N_mod_q[q_idx]) {
            return (candidate == q) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_candidate_prime_legendre_inline(uint64_t candidate) {
    int td = trial_division_legendre_inline(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_legendre_inline(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    precompute_N_mod(N);
    init_a_mod_inline(a_max);

    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_legendre_inline(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
        decrement_a_mod_inline();
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
    {"baseline",              find_solution_from_N},
    {"legendre (basic)",      find_solution_legendre},
    {"legendre (incr a%q)",   find_solution_legendre_fast},
    {"legendre (inline 3)",   find_solution_legendre_inline},
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

    printf("Legendre-Style Trial Division Test\n");
    printf("Iterations per scale: %s\n\n", fmt_num(count));

    /* Initialize sqrt tables */
    init_sqrt_tables();

    printf("Comparison:\n");
    printf("  Baseline: candidate %% q (64-bit modulo on large candidate)\n");
    printf("  Legendre: a %% q then square (can use incremental tracking)\n\n");

    /* Print header */
    printf("%-10s", "Scale");
    for (size_t s = 0; s < NUM_SOLVERS; s++) {
        printf("  %18s", SOLVERS[s].name);
    }
    printf("\n");
    for (size_t i = 0; i < 10 + NUM_SOLVERS * 20; i++) printf("-");
    printf("\n");

    /* Run benchmarks */
    for (size_t i = 0; i < NUM_SCALES; i++) {
        printf("%-10s", SCALES[i].label);

        double baseline_rate = 0;
        for (size_t s = 0; s < NUM_SOLVERS; s++) {
            double rate = run_benchmark(SOLVERS[s].func, SCALES[i].n_start, count);
            if (s == 0) baseline_rate = rate;

            double speedup = rate / baseline_rate;
            printf("  %10s (%4.2fx)", fmt_num((uint64_t)rate), speedup);
        }
        printf("\n");
    }

    return 0;
}
