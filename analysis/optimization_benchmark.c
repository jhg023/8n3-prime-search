/*
 * Optimization Benchmark: Compare different solution-finding strategies
 *
 * Tests:
 *   1. Baseline (current implementation)
 *   2. Root-class sieve (skip a where a² ≡ N mod q)
 *   3. Wheel/CRT (generate survivors directly)
 *   4. Combined optimizations
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/prime.h"
#include "../include/solve.h"

#include <inttypes.h>

/* ========================================================================== */
/* Timing utilities                                                           */
/* ========================================================================== */

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Tonelli-Shanks: Compute square root of n mod p (if it exists)              */
/* ========================================================================== */

/*
 * Returns r such that r² ≡ n (mod p), or 0 if no square root exists.
 * Assumes p is an odd prime.
 */
static uint64_t tonelli_shanks(uint64_t n, uint64_t p) {
    if (n == 0) return 0;
    if (p == 2) return n & 1;

    n %= p;
    if (n == 0) return 0;

    /* Check if n is a quadratic residue using Euler's criterion */
    /* n^((p-1)/2) ≡ 1 (mod p) iff n is a QR */
    if (powmod64(n, (p - 1) / 2, p) != 1) {
        return 0;  /* No square root exists */
    }

    /* Special case: p ≡ 3 (mod 4) */
    if ((p & 3) == 3) {
        return powmod64(n, (p + 1) / 4, p);
    }

    /* General Tonelli-Shanks algorithm */
    /* Write p - 1 = Q * 2^S with Q odd */
    uint64_t Q = p - 1;
    int S = 0;
    while ((Q & 1) == 0) {
        Q >>= 1;
        S++;
    }

    /* Find a quadratic non-residue z */
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

        /* Find the least i such that t^(2^i) ≡ 1 (mod p) */
        int i = 1;
        uint64_t temp = mulmod64(t, t, p);
        while (temp != 1) {
            temp = mulmod64(temp, temp, p);
            i++;
        }

        /* Update */
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
/* Strategy 1: Root-class sieve                                               */
/* ========================================================================== */

/*
 * Small primes for root-class sieve.
 * For each prime q, we'll compute r = sqrt(N) mod q (if it exists).
 * Then skip all a ≡ ±r (mod q).
 */
#define NUM_SIEVE_PRIMES 15
static const uint32_t SIEVE_PRIMES[NUM_SIEVE_PRIMES] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53
};

/*
 * Check if a should be skipped based on root-class sieve.
 * Returns true if the candidate p = (N - a²)/2 is guaranteed composite.
 */
static inline bool should_skip_a(uint64_t N, uint64_t a,
                                  const uint64_t* roots, const bool* has_root) {
    for (int i = 0; i < NUM_SIEVE_PRIMES; i++) {
        if (!has_root[i]) continue;

        uint32_t q = SIEVE_PRIMES[i];
        uint64_t a_mod = a % q;
        uint64_t r = roots[i];

        /* Skip if a ≡ r (mod q) or a ≡ -r (mod q) */
        /* But allow if candidate would equal q itself */
        if (a_mod == r || a_mod == (q - r) % q) {
            /* Check if candidate = q (the only case where divisible by q is OK) */
            uint64_t a_sq = a * a;
            if (a_sq <= N - 4) {
                uint64_t candidate = (N - a_sq) >> 1;
                if (candidate != q) {
                    return true;  /* Skip this a */
                }
            }
        }
    }
    return false;
}

static uint64_t find_solution_rootsieve(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute square roots of N mod each sieve prime */
    uint64_t roots[NUM_SIEVE_PRIMES];
    bool has_root[NUM_SIEVE_PRIMES];

    for (int i = 0; i < NUM_SIEVE_PRIMES; i++) {
        uint32_t q = SIEVE_PRIMES[i];
        uint64_t r = tonelli_shanks(N % q, q);
        roots[i] = r;
        has_root[i] = (r != 0 || N % q == 0);
    }

    uint64_t a = a_max;
    while (1) {
        /* Check if this a should be skipped */
        if (!should_skip_a(N, a, roots, has_root)) {
            uint64_t a_sq = a * a;

            if (a_sq <= N - 4) {
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
/* Strategy 2: Wheel/CRT - precompute valid residues                          */
/* ========================================================================== */

/*
 * Use a smaller wheel for practical implementation.
 * M = 3 * 5 * 7 = 105, so we track valid a mod 210 (since a must be odd)
 */
#define WHEEL_MOD 210  /* 2 * 3 * 5 * 7 */
#define WHEEL_PRIMES_COUNT 3
static const uint32_t WHEEL_PRIMES[WHEEL_PRIMES_COUNT] = {3, 5, 7};

/*
 * Precompute which residues of a (mod WHEEL_MOD) are valid for a given N.
 * A residue is invalid if for any wheel prime q, a² ≡ N (mod q).
 */
static int precompute_wheel(uint64_t N, uint8_t* valid_residues) {
    /* First, find roots of N mod each wheel prime */
    uint64_t roots[WHEEL_PRIMES_COUNT];
    bool has_root[WHEEL_PRIMES_COUNT];

    for (int i = 0; i < WHEEL_PRIMES_COUNT; i++) {
        uint32_t q = WHEEL_PRIMES[i];
        uint64_t r = tonelli_shanks(N % q, q);
        roots[i] = r;
        has_root[i] = (r != 0 || N % q == 0);
    }

    /* Mark valid residues (odd values only) */
    int count = 0;
    for (int a_mod = 1; a_mod < WHEEL_MOD; a_mod += 2) {
        bool valid = true;

        for (int i = 0; i < WHEEL_PRIMES_COUNT && valid; i++) {
            if (!has_root[i]) continue;

            uint32_t q = WHEEL_PRIMES[i];
            uint64_t r = roots[i];
            int am = a_mod % q;

            /* Invalid if a ≡ ±r (mod q) */
            if (am == (int)r || am == (int)((q - r) % q)) {
                valid = false;
            }
        }

        if (valid) {
            valid_residues[count++] = a_mod;
        }
    }

    return count;
}

static uint64_t find_solution_wheel(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute valid residues for this N */
    uint8_t valid_residues[WHEEL_MOD];
    int num_valid = precompute_wheel(N, valid_residues);

    if (num_valid == 0) {
        /* Extremely rare: no valid residues. Fall back to baseline. */
        return find_solution_baseline(n, p_out);
    }

    /* Iterate through valid residue classes in reverse order */
    /* Start from the largest a that fits in a valid residue class */
    uint64_t base = (a_max / WHEEL_MOD) * WHEEL_MOD;

    /* Find starting residue index */
    int start_idx = num_valid - 1;
    while (start_idx >= 0 && base + valid_residues[start_idx] > a_max) {
        start_idx--;
    }

    /* Main iteration */
    while (base > 0 || start_idx >= 0) {
        for (int i = start_idx; i >= 0; i--) {
            uint64_t a = base + valid_residues[i];

            if (a > a_max) continue;
            if (a < 1) continue;

            uint64_t a_sq = a * a;
            if (a_sq > N - 4) continue;

            uint64_t candidate = (N - a_sq) >> 1;

            /* Still need trial division for primes not in wheel */
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (base < WHEEL_MOD) break;
        base -= WHEEL_MOD;
        start_idx = num_valid - 1;
    }

    return 0;
}

/* ========================================================================== */
/* Strategy 3: Combined - root sieve with larger prime set + wheel            */
/* ========================================================================== */

/*
 * Combined approach:
 * 1. Use wheel mod 210 for fast iteration
 * 2. Apply root-class sieve for primes 11, 13, 17, 19, 23, 29, 31
 *    (primes 3, 5, 7 are already handled by the wheel)
 */

#define NUM_EXTRA_SIEVE_PRIMES 10
static const uint32_t EXTRA_SIEVE_PRIMES[NUM_EXTRA_SIEVE_PRIMES] = {
    11, 13, 17, 19, 23, 29, 31, 37, 41, 43
};

static inline bool should_skip_a_extra(uint64_t N, uint64_t a,
                                        const uint64_t* roots, const bool* has_root) {
    for (int i = 0; i < NUM_EXTRA_SIEVE_PRIMES; i++) {
        if (!has_root[i]) continue;

        uint32_t q = EXTRA_SIEVE_PRIMES[i];
        uint64_t a_mod = a % q;
        uint64_t r = roots[i];

        if (a_mod == r || a_mod == (q - r) % q) {
            uint64_t a_sq = a * a;
            if (a_sq <= N - 4) {
                uint64_t candidate = (N - a_sq) >> 1;
                if (candidate != q) {
                    return true;
                }
            }
        }
    }
    return false;
}

static uint64_t find_solution_combined(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    if ((a_max & 1) == 0) a_max--;

    /* Precompute wheel */
    uint8_t valid_residues[WHEEL_MOD];
    int num_valid = precompute_wheel(N, valid_residues);

    if (num_valid == 0) {
        return find_solution_baseline(n, p_out);
    }

    /* Precompute roots for extra sieve primes */
    uint64_t extra_roots[NUM_EXTRA_SIEVE_PRIMES];
    bool extra_has_root[NUM_EXTRA_SIEVE_PRIMES];

    for (int i = 0; i < NUM_EXTRA_SIEVE_PRIMES; i++) {
        uint32_t q = EXTRA_SIEVE_PRIMES[i];
        uint64_t r = tonelli_shanks(N % q, q);
        extra_roots[i] = r;
        extra_has_root[i] = (r != 0 || N % q == 0);
    }

    /* Iterate using wheel */
    uint64_t base = (a_max / WHEEL_MOD) * WHEEL_MOD;
    int start_idx = num_valid - 1;
    while (start_idx >= 0 && base + valid_residues[start_idx] > a_max) {
        start_idx--;
    }

    while (base > 0 || start_idx >= 0) {
        for (int i = start_idx; i >= 0; i--) {
            uint64_t a = base + valid_residues[i];

            if (a > a_max || a < 1) continue;

            /* Apply extra sieve */
            if (should_skip_a_extra(N, a, extra_roots, extra_has_root)) {
                continue;
            }

            uint64_t a_sq = a * a;
            if (a_sq > N - 4) continue;

            uint64_t candidate = (N - a_sq) >> 1;

            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (base < WHEEL_MOD) break;
        base -= WHEEL_MOD;
        start_idx = num_valid - 1;
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
    {"Root-sieve (15 primes)", find_solution_rootsieve},
    {"Wheel (mod 210)", find_solution_wheel},
    {"Combined (wheel + sieve)", find_solution_combined},
};

#define NUM_STRATEGIES (sizeof(strategies) / sizeof(strategies[0]))

static void benchmark_strategy(Strategy* s, uint64_t start_n, uint64_t count) {
    double t0 = get_time();
    uint64_t solutions = 0;
    uint64_t total_checks = 0;

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

    printf("  %-28s %10.0f n/sec  (%.3fs, %" PRIu64 " solutions)\n",
           s->name, rate, elapsed, solutions);
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
    printf("Optimization Benchmark: Root-class Sieve & Wheel\n");
    printf("=================================================\n");

    uint64_t count = 100000;  /* Default iteration count */

    if (argc > 1) {
        count = strtoull(argv[1], NULL, 10);
    }

    /* Test at multiple scales */
    uint64_t scales[] = {
        1000000000ULL,           /* 10^9 */
        1000000000000ULL,        /* 10^12 */
        1000000000000000ULL,     /* 10^15 */
        1000000000000000000ULL,  /* 10^18 */
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        run_benchmarks(scales[i], count);
    }

    /* Verify correctness with small range */
    printf("\nVerifying correctness (n = 1 to 10000)...\n");
    for (uint64_t n = 1; n <= 10000; n++) {
        uint64_t p_base, p_root, p_wheel, p_comb;
        uint64_t a_base = find_solution_baseline(n, &p_base);
        uint64_t a_root = find_solution_rootsieve(n, &p_root);
        uint64_t a_wheel = find_solution_wheel(n, &p_wheel);
        uint64_t a_comb = find_solution_combined(n, &p_comb);

        /* All should find a solution (no counterexamples known) */
        if (a_base == 0 || a_root == 0 || a_wheel == 0 || a_comb == 0) {
            printf("ERROR: Missing solution for n=%" PRIu64 "\n", n);
            printf("  baseline=%" PRIu64 ", rootsieve=%" PRIu64 ", wheel=%" PRIu64 ", combined=%" PRIu64 "\n",
                   a_base, a_root, a_wheel, a_comb);
            return 1;
        }

        /* Verify each found valid solutions */
        uint64_t N = 8 * n + 3;
        if (a_base * a_base + 2 * p_base != N ||
            a_root * a_root + 2 * p_root != N ||
            a_wheel * a_wheel + 2 * p_wheel != N ||
            a_comb * a_comb + 2 * p_comb != N) {
            printf("ERROR: Invalid solution for n=%" PRIu64 "\n", n);
            return 1;
        }
    }
    printf("All strategies produce valid solutions.\n");

    return 0;
}
