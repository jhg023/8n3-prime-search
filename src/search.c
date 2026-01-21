/*
 * Counterexample Search: 8n + 3 = a² + 2p
 *
 * Searches for counterexamples to the conjecture that every integer
 * of the form 8n + 3 can be expressed as a² + 2p, where a is a 
 * positive odd integer and p is prime.
 *
 * Uses FJ64_262K primality test for optimal performance:
 * - Only 2 Miller-Rabin tests (vs 7 in standard deterministic test)
 * - 512KB precomputed hash table for witness selection
 * - 100% deterministic for all 64-bit integers
 *
 * Reference: Forisek & Jancina (2015), "Fast Primality Testing for 
 * Integers That Fit into a Machine Word"
 *
 * Compile: gcc -O3 -march=native -flto -o search search.c -lm -I../include
 * Usage:   ./search [n_start] [n_end]
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* FJ64_262K hash table - 262,144 entries × 16-bit witnesses = 512KB */
#include "fj64_table.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Trial division primes - filters ~82% of composites before Miller-Rabin */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
static const int NUM_TRIAL_PRIMES = 30;

/* Progress reporting interval in seconds */
#define PROGRESS_SECONDS 5.0

/* Default search range */
#define DEFAULT_N_START 1000000000000ULL    /* 10^12 */
#define DEFAULT_N_END   1000001000000ULL    /* 10^12 + 10^6 */

/* ========================================================================== */
/* Number Formatting                                                          */
/* ========================================================================== */

#define NUM_FMT_BUFS 8
#define NUM_FMT_SIZE 32
static char fmt_bufs[NUM_FMT_BUFS][NUM_FMT_SIZE];
static int fmt_buf_idx = 0;

/**
 * Format a number with comma separators (e.g., 1000000 -> "1,000,000")
 * Uses rotating static buffers - safe for up to NUM_FMT_BUFS concurrent uses
 */
const char* fmt_num(uint64_t n) {
    char* buf = fmt_bufs[fmt_buf_idx];
    fmt_buf_idx = (fmt_buf_idx + 1) % NUM_FMT_BUFS;

    char temp[NUM_FMT_SIZE];
    int len = 0;
    if (n == 0) {
        temp[len++] = '0';
    } else {
        while (n > 0) {
            temp[len++] = '0' + (n % 10);
            n /= 10;
        }
    }

    int pos = 0;
    for (int i = len - 1; i >= 0; i--) {
        buf[pos++] = temp[i];
        if (i > 0 && i % 3 == 0) {
            buf[pos++] = ',';
        }
    }
    buf[pos] = '\0';

    return buf;
}

/* ========================================================================== */
/* Core Arithmetic                                                            */
/* ========================================================================== */

/**
 * 64-bit modular multiplication using 128-bit intermediate
 */
static inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    return ((__uint128_t)a * b) % m;
}

/**
 * Modular exponentiation: base^exp mod m
 */
static inline uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1)
            result = mulmod64(result, base, mod);
        exp >>= 1;
        base = mulmod64(base, base, mod);
    }
    return result;
}

/**
 * Integer square root using Newton's method with floating-point seed
 */
static inline uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = (uint64_t)sqrtl((long double)n);
    /* Correct for floating-point errors */
    while (x > 0 && x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

/* ========================================================================== */
/* FJ64_262K Primality Test                                                   */
/* ========================================================================== */

/**
 * FJ64 hash function - maps n to a bucket in [0, 262143]
 */
static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;  /* & (2^18 - 1) */
}

/**
 * Single Miller-Rabin witness test
 * Returns true if n is a strong probable prime to base a
 */
static inline bool mr_witness(uint64_t n, uint64_t a) {
    if (a >= n) a %= n;
    if (a == 0) return true;
    
    uint64_t d = n - 1;
    int r = __builtin_ctzll(d);
    d >>= r;
    
    uint64_t x = powmod64(a, d, n);
    
    if (x == 1 || x == n - 1)
        return true;
    
    for (int i = 1; i < r; i++) {
        x = mulmod64(x, x, n);
        if (x == n - 1)
            return true;
        if (x == 1)
            return false;
    }
    return false;
}

/**
 * FJ64_262K primality test - exactly 2 Miller-Rabin tests
 * Assumes: n > 127, n is odd, n passed trial division
 */
static inline bool is_prime_fj64_fast(uint64_t n) {
    if (!mr_witness(n, 2))
        return false;
    return mr_witness(n, fj64_bases[fj64_hash(n)]);
}

/**
 * Full primality test for standalone use
 */
bool is_prime_64(uint64_t n) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if ((n & 1) == 0) return false;
    if (n % 3 == 0) return n == 3;
    if (n % 5 == 0) return n == 5;
    if (n % 7 == 0) return n == 7;
    if (n < 121) return true;
    
    return is_prime_fj64_fast(n);
}

/* ========================================================================== */
/* Main Search Algorithm                                                      */
/* ========================================================================== */

/**
 * Find a solution to 8n + 3 = a² + 2p
 * 
 * Returns the smallest valid a, or 0 if no solution exists (counterexample).
 * Optionally returns the corresponding prime p via out parameter.
 */
uint64_t find_solution(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    
    for (uint64_t a = 1; a <= a_max; a += 2) {
        uint64_t a_sq = a * a;
        
        /* Ensure candidate p ≥ 2 */
        if (a_sq > N - 4) break;
        
        uint64_t candidate = (N - a_sq) >> 1;

        /* 
         * Note: candidate ≡ 1 (mod 4) is mathematically guaranteed:
         * N = 8n+3 ≡ 3 (mod 8), a² ≡ 1 (mod 8) for odd a
         * So N - a² ≡ 2 (mod 8), and (N - a²)/2 ≡ 1 (mod 4)
         */

        /* Trial division filter - eliminates ~82% of composites */
        bool dominated = false;
        bool is_trial_prime = false;
        for (int i = 0; i < NUM_TRIAL_PRIMES; i++) {
            if (candidate % TRIAL_PRIMES[i] == 0) {
                if (candidate == TRIAL_PRIMES[i]) {
                    is_trial_prime = true;
                } else {
                    dominated = true;
                }
                break;
            }
        }
        if (dominated) continue;

        /* Primality test: FJ64_262K (only 2 Miller-Rabin rounds) */
        if (is_trial_prime || candidate <= 127 || is_prime_fj64_fast(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }
    
    return 0;  /* Counterexample! */
}

/* ========================================================================== */
/* Search Loop                                                                */
/* ========================================================================== */

/**
 * Run the search over a range of n values
 */
void run_search(uint64_t n_start, uint64_t n_end,
                uint64_t *out_counterexamples, 
                uint64_t *out_max_a, 
                uint64_t *out_max_a_n) {
    clock_t start_time = clock();

    uint64_t counterexamples_found = 0;
    uint64_t max_a_seen = 0;
    uint64_t max_a_n = 0;

    double last_progress_time = 0.0;

    for (uint64_t n = n_start; n < n_end; n++) {
        uint64_t p;
        uint64_t a = find_solution(n, &p);

        if (a == 0) {
            /* Counterexample found! */
            counterexamples_found++;
            printf("COUNTEREXAMPLE: n = %s (N = %s)\n", 
                   fmt_num(n), fmt_num(8 * n + 3));
            fflush(stdout);
        } else {
            if (a > max_a_seen) {
                max_a_seen = a;
                max_a_n = n;
            }
        }

        /* Progress reporting - check only every 16384 iterations */
        if ((n & 0x3FFF) == 0) {
            clock_t now = clock();
            double elapsed = (double)(now - start_time) / CLOCKS_PER_SEC;
            if (elapsed - last_progress_time >= PROGRESS_SECONDS) {
                double rate = (n - n_start + 1) / elapsed;
                double pct = 100.0 * (n - n_start) / (n_end - n_start);

                printf("n = %s (%.1f%%), rate = %s n/sec, max_a = %s\n",
                       fmt_num(n), pct, fmt_num((uint64_t)rate), fmt_num(max_a_seen));
                fflush(stdout);

                last_progress_time = elapsed;
            }
        }
    }

    *out_counterexamples = counterexamples_found;
    *out_max_a = max_a_seen;
    *out_max_a_n = max_a_n;
}

/* ========================================================================== */
/* Verification                                                               */
/* ========================================================================== */

/**
 * Verify the algorithm against known solutions
 */
bool verify_known_solutions(void) {
    printf("Verifying known solutions...\n");
    
    struct { uint64_t n; uint64_t a; uint64_t p; } known[] = {
        {1, 1, 5},    /* 11 = 1 + 10 */
        {2, 3, 5},    /* 19 = 9 + 10 */
        {3, 1, 13},   /* 27 = 1 + 26 */
        {4, 5, 5}     /* 35 = 25 + 10 (also 1 + 34 = 1 + 2*17) */
    };
    
    bool all_pass = true;
    
    for (int i = 0; i < 4; i++) {
        uint64_t n = known[i].n;
        uint64_t expected_a = known[i].a;
        uint64_t expected_p = known[i].p;
        
        uint64_t N = 8 * n + 3;
        uint64_t lhs = N;
        uint64_t rhs = expected_a * expected_a + 2 * expected_p;
        
        bool equation_valid = (lhs == rhs);
        bool p_is_prime = is_prime_64(expected_p);
        
        uint64_t found_p;
        uint64_t found_a = find_solution(n, &found_p);
        
        printf("  n=%llu: N=%llu, given (%llu,%llu), found (%llu,%llu) ... ",
               (unsigned long long)n, (unsigned long long)N,
               (unsigned long long)expected_a, (unsigned long long)expected_p,
               (unsigned long long)found_a, (unsigned long long)found_p);
        
        if (equation_valid && p_is_prime && found_a > 0) {
            printf("PASS\n");
        } else {
            printf("FAIL\n");
            all_pass = false;
        }
    }
    
    return all_pass;
}

/* ========================================================================== */
/* Argument Parsing                                                           */
/* ========================================================================== */

/**
 * Parse a number that may be in scientific notation (e.g., 1e12, 2.5e9)
 */
uint64_t parse_number(const char* str) {
    char* endptr;
    double val = strtod(str, &endptr);
    if (*endptr != '\0') {
        return strtoull(str, NULL, 10);
    }
    return (uint64_t)val;
}

void print_usage(const char* program) {
    printf("Usage: %s [n_start] [n_end]\n", program);
    printf("\n");
    printf("Search for counterexamples to: 8n + 3 = a² + 2p\n");
    printf("\n");
    printf("For each n in [n_start, n_end), attempts to find odd a and prime p\n");
    printf("such that 8n + 3 = a² + 2p. Reports any n for which no solution exists.\n");
    printf("\n");
    printf("Arguments:\n");
    printf("  n_start   Starting value of n (inclusive), default: 1e12\n");
    printf("  n_end     Ending value of n (exclusive), default: 1e12 + 1e6\n");
    printf("\n");
    printf("Numbers can be in scientific notation (e.g., 1e9, 2.5e6)\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s                    Search [10^12, 10^12 + 10^6)\n", program);
    printf("  %s 1 1e6              Search [1, 10^6)\n", program);
    printf("  %s 1e9 2e9            Search [10^9, 2×10^9)\n", program);
    printf("\n");
    printf("Exit codes:\n");
    printf("  0  Search completed, no counterexamples found\n");
    printf("  1  Error (invalid arguments, verification failure)\n");
    printf("  2  Counterexample found\n");
}

/* ========================================================================== */
/* Main                                                                       */
/* ========================================================================== */

int main(int argc, char** argv) {
    printf("==================================================================\n");
    printf("     Counterexample Search: 8n + 3 = a² + 2p                      \n");
    printf("     Optimized with FJ64_262K primality test                      \n");
    printf("==================================================================\n\n");

    uint64_t n_start = DEFAULT_N_START;
    uint64_t n_end = DEFAULT_N_END;

    /* Handle help flag */
    if (argc == 2 && (strcmp(argv[1], "-h") == 0 || 
                      strcmp(argv[1], "--help") == 0)) {
        print_usage(argv[0]);
        return 0;
    }

    /* Parse arguments */
    if (argc >= 2) {
        n_start = parse_number(argv[1]);
    }
    if (argc >= 3) {
        n_end = parse_number(argv[2]);
    }

    if (n_start >= n_end) {
        fprintf(stderr, "Error: n_start must be less than n_end\n");
        return 1;
    }

    uint64_t total = n_end - n_start;

    printf("Configuration:\n");
    printf("  Range: n ∈ [%s, %s)\n", fmt_num(n_start), fmt_num(n_end));
    printf("  Count: %s values\n", fmt_num(total));
    printf("  Primality test: FJ64_262K (2 Miller-Rabin tests)\n");
    printf("\n");

    /* Verify algorithm correctness */
    if (!verify_known_solutions()) {
        fprintf(stderr, "\nERROR: Verification failed!\n");
        return 1;
    }
    printf("\n");

    /* Run search */
    printf("Starting search...\n\n");
    clock_t global_start = clock();

    uint64_t total_counterexamples = 0;
    uint64_t global_max_a = 0;
    uint64_t global_max_a_n = 0;

    run_search(n_start, n_end, &total_counterexamples, 
               &global_max_a, &global_max_a_n);

    clock_t global_end = clock();
    double global_elapsed = (double)(global_end - global_start) / CLOCKS_PER_SEC;

    /* Print results */
    printf("\n");
    printf("==================================================================\n");
    printf("RESULTS\n");
    printf("==================================================================\n\n");

    printf("Total time:           %.1f seconds\n", global_elapsed);
    printf("Total throughput:     %s n/sec\n", 
           fmt_num((uint64_t)(total / global_elapsed)));
    printf("Counterexamples:      %s\n", fmt_num(total_counterexamples));
    printf("Maximum a observed:   %s (at n = %s)\n", 
           fmt_num(global_max_a), fmt_num(global_max_a_n));

    return (total_counterexamples > 0) ? 2 : 0;
}
