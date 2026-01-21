/*
 * Analyze prime candidate sizes for different a values
 * Shows what size primes are tested when iterating large-to-small
 */
#include <stdio.h>
#include <stdint.h>
#include <math.h>

static inline uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = (uint64_t)sqrt((double)n);
    while (x * x > n) x--;
    while ((x + 1) * (x + 1) <= n) x++;
    return x;
}

int main(void) {
    printf("=================================================================\n");
    printf("Prime candidate sizes for 8n+3 = a² + 2p (large-to-small iteration)\n");
    printf("=================================================================\n\n");

    uint64_t test_n[] = {
        1000000ULL,           // 10^6
        1000000000ULL,        // 10^9
        1000000000000ULL,     // 10^12
        1000000000000000ULL   // 10^15
    };
    const char* labels[] = {"10^6", "10^9", "10^12", "10^15"};

    for (int t = 0; t < 4; t++) {
        uint64_t n = test_n[t];
        uint64_t N = 8 * n + 3;
        uint64_t a_max = isqrt64(N);
        if ((a_max & 1) == 0) a_max--;

        printf("n = %s (N = 8n+3 = %llu)\n", labels[t], (unsigned long long)N);
        printf("a_max = %llu\n", (unsigned long long)a_max);
        printf("--------------------------------------------------\n");
        printf("  %-12s  %-20s  %-8s\n", "a value", "p = (N-a²)/2", "bits");
        printf("--------------------------------------------------\n");

        // Show first 10 a values (largest a, smallest p)
        uint64_t count = 0;
        for (uint64_t a = a_max; a >= 1 && count < 10; a -= 2) {
            uint64_t a_sq = a * a;
            if (a_sq > N - 4) continue;
            uint64_t p = (N - a_sq) >> 1;
            int bits = 64 - __builtin_clzll(p);
            printf("  %-12llu  %-20llu  %d bits\n",
                   (unsigned long long)a, (unsigned long long)p, bits);
            count++;
        }

        printf("  ...\n");

        // Show last 5 a values (smallest a, largest p)
        uint64_t small_a[5];
        uint64_t small_p[5];
        int idx = 0;
        for (uint64_t a = 1; a <= a_max && idx < 5; a += 2) {
            uint64_t a_sq = a * a;
            if (a_sq > N - 4) break;
            small_a[idx] = a;
            small_p[idx] = (N - a_sq) >> 1;
            idx++;
        }
        for (int i = 0; i < idx; i++) {
            int bits = 64 - __builtin_clzll(small_p[i]);
            printf("  %-12llu  %-20llu  %d bits\n",
                   (unsigned long long)small_a[i], (unsigned long long)small_p[i], bits);
        }

        // Calculate statistics
        uint64_t p_at_max_a = (N - a_max * a_max) >> 1;
        uint64_t p_at_a1 = (N - 1) >> 1;

        printf("\nSummary:\n");
        printf("  Largest a tested:  a=%llu → p=%llu (%d bits)\n",
               (unsigned long long)a_max, (unsigned long long)p_at_max_a,
               64 - __builtin_clzll(p_at_max_a));
        printf("  Smallest a tested: a=1 → p=%llu (%d bits)\n",
               (unsigned long long)p_at_a1, 64 - __builtin_clzll(p_at_a1));

        // Count how many a values have sievable p
        uint64_t sievable_1m = 0, sievable_10m = 0, sievable_100m = 0;
        uint64_t total_a = 0;
        for (uint64_t a = a_max; a >= 1; a -= 2) {
            uint64_t a_sq = a * a;
            if (a_sq > N - 4) continue;
            total_a++;
            uint64_t p = (N - a_sq) >> 1;
            if (p < 1000000ULL) sievable_1m++;
            if (p < 10000000ULL) sievable_10m++;
            if (p < 100000000ULL) sievable_100m++;
        }

        printf("\nSieve potential (iterating large-to-small):\n");
        printf("  p < 1M:   %llu/%llu a values (%.1f%%)\n",
               (unsigned long long)sievable_1m, (unsigned long long)total_a,
               100.0 * sievable_1m / total_a);
        printf("  p < 10M:  %llu/%llu a values (%.1f%%)\n",
               (unsigned long long)sievable_10m, (unsigned long long)total_a,
               100.0 * sievable_10m / total_a);
        printf("  p < 100M: %llu/%llu a values (%.1f%%)\n",
               (unsigned long long)sievable_100m, (unsigned long long)total_a,
               100.0 * sievable_100m / total_a);

        printf("\n\n");
    }

    return 0;
}
