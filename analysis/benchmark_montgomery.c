/*
 * Benchmark: Montgomery vs __uint128_t for Miller-Rabin
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "arith.h"
#include "arith_montgomery.h"
#include "fj64_table.h"

/* FJ64 hash function */
static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

/* Standard Miller-Rabin witness (from prime.h) */
static inline bool mr_witness_standard(uint64_t n, uint64_t a) {
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

/* FJ64 test using standard arithmetic */
static inline bool is_prime_fj64_standard(uint64_t n) {
    if (!mr_witness_standard(n, 2))
        return false;
    return mr_witness_standard(n, fj64_bases[fj64_hash(n)]);
}

/* FJ64 test using Montgomery arithmetic */
static inline bool is_prime_fj64_montgomery(uint64_t n) {
    if (!mr_witness_montgomery(n, 2))
        return false;
    return mr_witness_montgomery(n, fj64_bases[fj64_hash(n)]);
}

/* Generate prime candidates (simple sieve for small values) */
static inline bool is_candidate_for_bench(uint64_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;
    if (n % 3 == 0) return n == 3;
    if (n % 5 == 0) return n == 5;
    if (n % 7 == 0) return n == 7;
    return n > 127;
}

int main(void) {
    printf("Montgomery vs Standard Benchmark\n");
    printf("================================\n\n");

    /* Test at different scales */
    uint64_t starts[] = {
        1000000000ULL,            /* 10^9 */
        1000000000000ULL,         /* 10^12 */
        1000000000000000ULL,      /* 10^15 */
        2000000000000000000ULL    /* 2e18 */
    };
    const char* labels[] = {"10^9", "10^12", "10^15", "2e18"};

    for (int t = 0; t < 4; t++) {
        uint64_t start = starts[t];
        uint64_t count = 100000;

        /* Generate candidates that would pass trial division */
        uint64_t candidates[10000];
        int num_candidates = 0;

        for (uint64_t n = start; num_candidates < 10000; n++) {
            if (is_candidate_for_bench(n)) {
                candidates[num_candidates++] = n;
            }
        }

        printf("Testing at ~%s with %d candidates\n", labels[t], num_candidates);

        /* Benchmark standard */
        volatile int prime_count_std = 0;
        clock_t t0 = clock();
        for (int i = 0; i < num_candidates; i++) {
            if (is_prime_fj64_standard(candidates[i])) {
                prime_count_std++;
            }
        }
        clock_t t1 = clock();
        double time_std = (double)(t1 - t0) / CLOCKS_PER_SEC;

        /* Benchmark Montgomery */
        volatile int prime_count_mont = 0;
        clock_t t2 = clock();
        for (int i = 0; i < num_candidates; i++) {
            if (is_prime_fj64_montgomery(candidates[i])) {
                prime_count_mont++;
            }
        }
        clock_t t3 = clock();
        double time_mont = (double)(t3 - t2) / CLOCKS_PER_SEC;

        /* Verify results match */
        if (prime_count_std != prime_count_mont) {
            printf("  ERROR: Results mismatch! std=%d, mont=%d\n",
                   prime_count_std, prime_count_mont);
        }

        printf("  Standard:   %.4f sec (%.0f tests/sec)\n",
               time_std, num_candidates / time_std);
        printf("  Montgomery: %.4f sec (%.0f tests/sec)\n",
               time_mont, num_candidates / time_mont);
        printf("  Speedup:    %.2fx\n", time_std / time_mont);
        printf("  Primes found: %d\n\n", prime_count_std);
    }

    /* Verify correctness on known primes and composites */
    printf("Correctness verification:\n");

    uint64_t test_primes[] = {
        2, 3, 5, 7, 11, 13, 127, 131,
        1000000007, 1000000009,
        18446744073709551557ULL  /* Largest 64-bit prime */
    };
    uint64_t test_composites[] = {
        4, 9, 15, 121, 1000000008,
        18446744073709551615ULL  /* 2^64 - 1 = 3 * 5 * 17 * ... */
    };

    bool all_pass = true;

    for (int i = 0; i < 11; i++) {
        uint64_t n = test_primes[i];
        if (n < 128) continue;  /* Skip small primes handled specially */
        bool std = is_prime_fj64_standard(n);
        bool mont = is_prime_fj64_montgomery(n);
        if (std != mont || !std) {
            printf("  FAIL: %llu should be prime (std=%d, mont=%d)\n",
                   (unsigned long long)n, std, mont);
            all_pass = false;
        }
    }

    for (int i = 0; i < 6; i++) {
        uint64_t n = test_composites[i];
        if (n < 128) continue;
        if ((n & 1) == 0) continue;  /* Skip even (not tested) */
        bool std = is_prime_fj64_standard(n);
        bool mont = is_prime_fj64_montgomery(n);
        if (std != mont || std) {
            printf("  FAIL: %llu should be composite (std=%d, mont=%d)\n",
                   (unsigned long long)n, std, mont);
            all_pass = false;
        }
    }

    if (all_pass) {
        printf("  All correctness tests passed!\n");
    }

    return 0;
}
