/*
 * Test and benchmark prime_sieve_fast.h
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

/* Include the fast sieve */
#include "../include/prime_sieve_fast.h"

/* Reference primality test for verification */
static bool is_prime_ref(uint64_t n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if ((n & 1) == 0) return false;
    for (uint64_t i = 3; i * i <= n; i += 2) {
        if (n % i == 0) return false;
    }
    return true;
}

static double get_time_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1000000.0;
}

int main(int argc, char *argv[]) {
    uint64_t threshold = 1000000000ULL;  /* 10^9 default */
    if (argc > 1) {
        threshold = strtoull(argv[1], NULL, 10);
    }

    printf("Fast Sieve Test\n");
    printf("===============\n");
    printf("Threshold: %llu (%.2e)\n\n", (unsigned long long)threshold, (double)threshold);

#ifdef _OPENMP
    printf("OpenMP threads: %d\n\n", omp_get_max_threads());
#else
    printf("OpenMP: disabled\n\n");
#endif

    /* Create sieve */
    printf("Creating sieve...\n");
    double t0 = get_time_ms();
    PrimeSieve *sieve = sieve_create(threshold);
    double t1 = get_time_ms();

    if (!sieve) {
        fprintf(stderr, "Failed to create sieve\n");
        return 1;
    }

    printf("  Creation time: %.1f ms\n", t1 - t0);
    printf("  Memory: %s\n\n", sieve_memory_str(sieve));

    /* Verify correctness for small primes */
    printf("Verifying correctness...\n");
    int errors = 0;

    /* Test all numbers up to 10000 */
    for (uint64_t n = 0; n <= 10000 && n <= threshold; n++) {
        bool expected = is_prime_ref(n);
        bool got = sieve_is_prime(sieve, n);
        if (got != expected) {
            printf("  ERROR: n=%llu expected=%d got=%d\n",
                   (unsigned long long)n, expected, got);
            errors++;
            if (errors > 10) {
                printf("  ... too many errors, stopping\n");
                break;
            }
        }
    }

    /* Test some known primes */
    uint64_t test_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                              97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
                              1009, 10007, 100003, 1000003, 10000019, 100000007};
    for (int i = 0; i < (int)(sizeof(test_primes)/sizeof(test_primes[0])); i++) {
        uint64_t p = test_primes[i];
        if (p > threshold) continue;
        if (!sieve_is_prime(sieve, p)) {
            printf("  ERROR: %llu should be prime\n", (unsigned long long)p);
            errors++;
        }
    }

    /* Test some known composites */
    uint64_t test_composites[] = {4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20, 21, 22,
                                   25, 26, 27, 28, 33, 35, 49, 77, 91, 121, 143,
                                   1001, 10001, 100001, 1000001};
    for (int i = 0; i < (int)(sizeof(test_composites)/sizeof(test_composites[0])); i++) {
        uint64_t c = test_composites[i];
        if (c > threshold) continue;
        if (sieve_is_prime(sieve, c)) {
            printf("  ERROR: %llu should be composite\n", (unsigned long long)c);
            errors++;
        }
    }

    if (errors == 0) {
        printf("  All tests passed ✓\n\n");
    } else {
        printf("  %d errors found ✗\n\n", errors);
    }

    /* Count primes */
    printf("Counting primes...\n");
    t0 = get_time_ms();
    uint64_t count = sieve_prime_count(sieve);
    t1 = get_time_ms();
    printf("  Prime count: %llu\n", (unsigned long long)count);
    printf("  Count time: %.1f ms\n\n", t1 - t0);

    /* Expected counts for verification */
    /* pi(10^8) = 5,761,455 */
    /* pi(10^9) = 50,847,534 */
    /* pi(10^10) = 455,052,511 */
    if (threshold == 100000000ULL && count != 5761455) {
        printf("  WARNING: Expected 5761455 for 10^8\n");
    } else if (threshold == 1000000000ULL && count != 50847534) {
        printf("  WARNING: Expected 50847534 for 10^9\n");
    } else if (threshold == 10000000000ULL && count != 455052511) {
        printf("  WARNING: Expected 455052511 for 10^10\n");
    }

    /* Benchmark lookup speed */
    printf("Benchmarking lookup speed...\n");
    uint64_t lookups = 10000000;
    uint64_t found = 0;
    t0 = get_time_ms();
    for (uint64_t i = 0; i < lookups; i++) {
        uint64_t n = (i * 1000003) % threshold + 1;
        if (sieve_is_prime(sieve, n)) found++;
    }
    t1 = get_time_ms();
    printf("  %llu lookups in %.1f ms (%.1f M lookups/sec)\n",
           (unsigned long long)lookups, t1 - t0,
           lookups / (t1 - t0) / 1000.0);
    printf("  Found %llu primes (sanity check)\n\n", (unsigned long long)found);

    sieve_destroy(sieve);

    printf("===============\n");
    printf("Done.\n");

    return errors > 0 ? 1 : 0;
}
