/*
 * Test file for Legendre symbol table optimization
 *
 * The optimization idea: For small primes q, precompute sqrt[q][r] tables
 * that give a square root of r mod q (or NONE if r is not a QR).
 *
 * This allows instant Legendre symbol computation: (N/q) = -1 iff sqrt[q][N%q] == NONE
 *
 * Intended use case: "hard mode" when simple trial division doesn't find
 * a prime quickly - we can use Legendre symbols to skip impossible candidates.
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
/* Legendre Tables                                                            */
/* ========================================================================== */

/*
 * Small primes for Legendre testing
 * Using primes up to 661 as suggested (first 120 odd primes)
 */
static const uint16_t LEGENDRE_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157,
    163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239,
    241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331,
    337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421,
    431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509,
    521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613,
    617, 619, 631, 641, 643, 647, 653, 659, 661
};
#define NUM_LEGENDRE_PRIMES 120
#define MAX_LEGENDRE_PRIME 661

#define SQRT_NONE 0xFF

/*
 * sqrt_table[q_idx][r] = sqrt(r) mod LEGENDRE_PRIMES[q_idx], or SQRT_NONE if not QR
 * For q=661 (largest), we need 661 bytes per prime = ~79KB total
 */
static uint8_t sqrt_table[NUM_LEGENDRE_PRIMES][MAX_LEGENDRE_PRIME + 1];

/*
 * is_qr_table[q_idx][r] = 1 if r is a QR mod LEGENDRE_PRIMES[q_idx], 0 otherwise
 * More compact than sqrt_table if we only need QR testing
 */
static uint8_t is_qr_table[NUM_LEGENDRE_PRIMES][MAX_LEGENDRE_PRIME + 1];

/*
 * Initialize the square root and QR tables
 */
static void init_legendre_tables(void) {
    memset(sqrt_table, SQRT_NONE, sizeof(sqrt_table));
    memset(is_qr_table, 0, sizeof(is_qr_table));

    for (int q_idx = 0; q_idx < NUM_LEGENDRE_PRIMES; q_idx++) {
        uint16_t q = LEGENDRE_PRIMES[q_idx];

        /* 0 is always a QR with sqrt = 0 */
        sqrt_table[q_idx][0] = 0;
        is_qr_table[q_idx][0] = 1;

        /* Compute squares: x^2 mod q for x in [1, (q-1)/2] */
        for (uint16_t x = 1; x <= (q - 1) / 2; x++) {
            uint16_t sq = (x * x) % q;
            sqrt_table[q_idx][sq] = (uint8_t)x;
            is_qr_table[q_idx][sq] = 1;
        }
    }
}

/*
 * Check if r is a quadratic residue mod q using table lookup
 */
static inline bool is_qr_fast(int q_idx, uint64_t r) {
    uint16_t q = LEGENDRE_PRIMES[q_idx];
    return is_qr_table[q_idx][r % q];
}

/*
 * Get a square root of r mod q, or SQRT_NONE if not a QR
 */
static inline uint8_t sqrt_mod_fast(int q_idx, uint64_t r) {
    uint16_t q = LEGENDRE_PRIMES[q_idx];
    return sqrt_table[q_idx][r % q];
}

/* ========================================================================== */
/* Alternative find_solution with Legendre filtering                          */
/* ========================================================================== */

/*
 * Bitset sieve approach: For K candidates, compute a bitmask of which ones
 * are "bad" (divisible by some small q), then only run primality on survivors.
 */
#define SIEVE_K 32

/*
 * Find solution using Legendre-based filtering in "hard mode"
 *
 * Strategy:
 * 1. First try simple trial division (same as current)
 * 2. If no quick solution found in first few tries, switch to "hard mode"
 * 3. In hard mode, use Legendre symbols to filter candidates
 */
static inline uint64_t find_solution_legendre(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    /* Try normal approach first for QUICK_TRIES iterations */
    #define QUICK_TRIES 8

    for (int tries = 0; tries < QUICK_TRIES && a >= 1; tries++) {
        if (candidate >= 2) {
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    /*
     * "Hard mode": We didn't find a quick solution.
     *
     * Now use Legendre filtering. For each remaining candidate p = (N - a²)/2,
     * we want p to be prime. But if p has (p/q) = -1 for some small q,
     * and p % q != 0, then p is definitely composite IF q² > p (since
     * then q cannot divide p evenly).
     *
     * Actually, the simpler approach: just use extended trial division
     * with more primes, which is what the bitset sieve effectively does.
     *
     * The "table-driven Legendre" optimization is more useful for checking
     * whether N itself has a representation, not for filtering individual
     * candidates. Let me implement what the image actually suggests:
     *
     * Use the sqrt tables to find which 'a' values could work mod each q.
     */

    /* Continue with normal iteration for remaining candidates */
    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    return 0;  /* Counterexample */
}

/*
 * Bitset sieve approach from the image
 *
 * For K candidates starting from a_max, compute which are "bad" (have p
 * divisible by some small q), then only test survivors.
 */
static inline uint64_t find_solution_bitset_sieve(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    /* Batch of K candidates */
    uint64_t a_values[SIEVE_K];
    uint64_t candidates[SIEVE_K];
    uint32_t bad_mask = 0;

    /* Initialize the first K candidates */
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    int k = 0;
    while (k < SIEVE_K && a >= 1) {
        if (candidate >= 2) {
            a_values[k] = a;
            candidates[k] = candidate;
            k++;
        }
        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    if (k == 0) return 0;

    /*
     * Sieve: For each small prime q, mark candidates where candidate % q == 0
     * (and candidate != q, meaning composite)
     */
    for (int i = 0; i < NUM_TRIAL_PRIMES && bad_mask != (1U << k) - 1; i++) {
        uint32_t q = TRIAL_PRIMES[i];
        for (int j = 0; j < k; j++) {
            if (!(bad_mask & (1U << j))) {
                if (candidates[j] % q == 0 && candidates[j] != q) {
                    bad_mask |= (1U << j);
                }
            }
        }
    }

    /* Test survivors with full primality */
    for (int j = 0; j < k; j++) {
        if (!(bad_mask & (1U << j))) {
            uint64_t cand = candidates[j];
            /* Check if it's actually prime */
            if (cand <= 127) {
                /* Small primes already passed trial division */
                if (p_out) *p_out = cand;
                return a_values[j];
            }
            if (is_prime_fj64_fast(cand)) {
                if (p_out) *p_out = cand;
                return a_values[j];
            }
        }
    }

    /* Continue with remaining candidates using normal approach */
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

/*
 * Extended trial division with more primes (up to 661)
 * This is a simpler optimization than the Legendre approach
 */
static inline int trial_division_extended(uint64_t candidate) {
    /* Check against all 120 small primes */
    for (int i = 0; i < NUM_LEGENDRE_PRIMES; i++) {
        uint16_t q = LEGENDRE_PRIMES[i];
        if (candidate % q == 0) {
            return (candidate == q) ? 1 : 0;  /* 1 = is prime q, 0 = composite */
        }
    }
    return 2;  /* Needs Miller-Rabin */
}

static inline bool is_candidate_prime_extended(uint64_t candidate) {
    int td = trial_division_extended(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    return is_prime_fj64_fast(candidate);
}

static inline uint64_t find_solution_extended_td(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_extended(candidate)) {
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
    {"baseline (current)",     find_solution_from_N},
    {"extended TD (120 primes)", find_solution_extended_td},
    {"bitset sieve (K=32)",    find_solution_bitset_sieve},
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

    printf("Legendre/Extended Trial Division Optimization Test\n");
    printf("Iterations per scale: %s\n\n", fmt_num(count));

    /* Initialize Legendre tables */
    printf("Initializing Legendre tables...\n");
    init_legendre_tables();
    printf("Tables initialized (%.1f KB)\n\n",
           (sizeof(sqrt_table) + sizeof(is_qr_table)) / 1024.0);

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
            double rate = run_benchmark(SOLVERS[s].func, SCALES[i].n_start, count);
            if (s == 0) baseline_rate = rate;

            double speedup = rate / baseline_rate;
            printf("  %12s (%4.2fx)", fmt_num((uint64_t)rate), speedup);
        }
        printf("\n");
    }

    printf("\nNote: Extended TD uses 120 primes (up to 661) vs baseline 30 primes (up to 127)\n");
    printf("      Bitset sieve batches %d candidates and sieves with %d primes\n",
           SIEVE_K, NUM_TRIAL_PRIMES);

    return 0;
}
