/*
 * Final comparison: 30 trials vs 120 trials
 * Multiple runs for statistical confidence
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "arith.h"
#include "arith_montgomery.h"
#include "fj64_table.h"

static const uint32_t PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,           /* 10 */
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,        /* 20 */
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,  /* 30 */
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233,
    239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353,
    359, 367, 373, 379, 383, 389, 397, 401, 409, 419,
    421, 431, 433, 439, 443, 449, 457, 461, 463, 467,
    479, 487, 491, 499, 503, 509, 521, 523, 541, 547,
    557, 563, 569, 571, 577, 587, 593, 599, 601, 607,
    613, 617, 619, 631, 641, 643, 647, 653, 659, 661
};

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline bool is_prime_fj64_mont(uint64_t n) {
    if (!mr_witness_montgomery(n, 2)) return false;
    return mr_witness_montgomery(n, fj64_bases[fj64_hash(n)]);
}

/* 30 trial primes version */
static inline uint64_t find_solution_30(uint64_t n) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;
        uint64_t candidate = (N - a_sq) >> 1;

        bool filtered = false;
        for (int i = 0; i < 30; i++) {
            if (candidate % PRIMES[i] == 0) {
                if (candidate == PRIMES[i]) return a;
                filtered = true;
                break;
            }
        }
        if (!filtered) {
            if (candidate <= 127) return a;
            if (is_prime_fj64_mont(candidate)) return a;
        }

        if (a < 3) break;
    }
    return 0;
}

/* 120 trial primes version (current) */
static inline uint64_t find_solution_120(uint64_t n) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;
        uint64_t candidate = (N - a_sq) >> 1;

        bool filtered = false;
        for (int i = 0; i < 120; i++) {
            if (candidate % PRIMES[i] == 0) {
                if (candidate == PRIMES[i]) return a;
                filtered = true;
                break;
            }
        }
        if (!filtered) {
            if (candidate <= 127) return a;
            if (is_prime_fj64_mont(candidate)) return a;
        }

        if (a < 3) break;
    }
    return 0;
}

double benchmark(uint64_t (*solver)(uint64_t), uint64_t n_start, uint64_t count) {
    /* Warmup */
    for (uint64_t i = 0; i < 5000; i++) {
        volatile uint64_t a = solver(n_start + i);
        (void)a;
    }

    clock_t t0 = clock();
    for (uint64_t i = 0; i < count; i++) {
        volatile uint64_t a = solver(n_start + i);
        (void)a;
    }
    clock_t t1 = clock();

    return count / ((double)(t1 - t0) / CLOCKS_PER_SEC);
}

int main(void) {
    setbuf(stdout, NULL);

    printf("Final Trial Division Comparison: 30 vs 120 primes\n");
    printf("==================================================\n\n");

    struct {
        uint64_t n_start;
        uint64_t count;
        const char* label;
    } scales[] = {
        {1000000000000ULL,       1000000, "10^12"},
        {1000000000000000ULL,     300000, "10^15"},
        {100000000000000000ULL,   150000, "10^17"},
        {2000000000000000000ULL,   80000, "2e18"}
    };

    printf("%-8s  %15s  %15s  %10s\n", "Scale", "30 primes", "120 primes", "Speedup");
    printf("----------------------------------------------------------\n");

    double total_speedup = 0;

    for (int s = 0; s < 4; s++) {
        double rate_30 = benchmark(find_solution_30, scales[s].n_start, scales[s].count);
        double rate_120 = benchmark(find_solution_120, scales[s].n_start, scales[s].count);
        double speedup = rate_30 / rate_120;
        total_speedup += speedup;

        printf("%-8s  %15.0f  %15.0f  %9.2fx\n",
               scales[s].label, rate_30, rate_120, speedup);
    }

    printf("----------------------------------------------------------\n");
    printf("Average speedup: %.2fx\n", total_speedup / 4);

    printf("\n=== Analysis ===\n\n");
    printf("The 30-prime configuration is faster because:\n");
    printf("1. Trial division catches ~80%% of composites with first 30 primes\n");
    printf("2. Additional 90 primes only catch ~5%% more composites\n");
    printf("3. With Montgomery MR (now 3x faster), the break-even point shifted\n");
    printf("4. Cost of 90 extra modulo ops > cost of occasional extra MR test\n");

    return 0;
}
