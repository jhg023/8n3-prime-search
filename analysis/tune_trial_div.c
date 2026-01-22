/*
 * Tune trial division prime count for different scales
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <inttypes.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/prime.h"

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* Extended prime list */
static const uint32_t ALL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,           /* 10 */
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,         /* 20 */
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127,   /* 30 */
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, /* 40 */
    181, 191, 193, 197, 199, 211, 223, 227, 229, 233  /* 50 */
};

/* Parameterized trial division */
static inline int trial_div_n(uint64_t candidate, int num_primes) {
    for (int i = 0; i < num_primes; i++) {
        if (candidate % ALL_PRIMES[i] == 0) {
            return (candidate == ALL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool is_prime_n(uint64_t candidate, int num_primes) {
    int td = trial_div_n(candidate, num_primes);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= ALL_PRIMES[num_primes - 1]) return true;
    return is_prime_fj64_fast(candidate);
}

static uint64_t find_solution_n(uint64_t n, uint64_t* p_out, int num_primes) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;

        uint64_t candidate = (N - a_sq) >> 1;
        if (is_prime_n(candidate, num_primes)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }
    return 0;
}

static void benchmark(uint64_t start_n, uint64_t count, int num_primes) {
    double t0 = get_time();

    volatile uint64_t sum = 0;  /* Prevent optimizer from eliminating loop */
    for (uint64_t n = start_n; n < start_n + count; n++) {
        uint64_t p;
        uint64_t a = find_solution_n(n, &p, num_primes);
        sum += a;
    }

    double elapsed = get_time() - t0;
    if (elapsed < 0.001) elapsed = 0.001;  /* Avoid division by zero */
    printf("    %2d primes: %10.0f n/sec  (%.3fs, sum=%" PRIu64 ")\n",
           num_primes, (double)count / elapsed, elapsed, sum);
}

int main(int argc, char** argv) {
    printf("Trial Division Prime Count Tuning\n");
    printf("==================================\n\n");

    uint64_t count = 100000;
    if (argc > 1) count = strtoull(argv[1], NULL, 10);

    int prime_counts[] = {10, 15, 20, 25, 30, 35, 40, 50};
    int num_counts = sizeof(prime_counts) / sizeof(prime_counts[0]);

    uint64_t scales[] = {
        1000000000ULL,
        1000000000000ULL,
        1000000000000000ULL,
        1000000000000000000ULL,
    };

    for (size_t s = 0; s < sizeof(scales) / sizeof(scales[0]); s++) {
        printf("Scale: %.0e\n", (double)scales[s]);
        for (int i = 0; i < num_counts; i++) {
            benchmark(scales[s], count, prime_counts[i]);
        }
        printf("\n");
    }

    return 0;
}
