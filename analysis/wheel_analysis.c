/*
 * Analyze which a values produce candidates divisible by small primes
 * Goal: Find patterns to skip a values that can never yield primes
 */
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

int main(void) {
    printf("Wheel Analysis for 8n + 3 = a² + 2p\n");
    printf("====================================\n\n");

    /*
     * For N = 8n + 3, candidate p = (N - a²) / 2
     *
     * We want to find when candidate is divisible by small primes.
     * This depends on (n mod m) and (a mod 2m) for prime m.
     */

    printf("=== Divisibility by 3 ===\n\n");
    printf("N = 8n + 3, candidate = (N - a²) / 2\n");
    printf("candidate mod 3 = ((8n + 3 - a²) / 2) mod 3\n\n");

    printf("Truth table: candidate mod 3 for each (n mod 3, a mod 6):\n");
    printf("  n mod 3 | a≡1(6) | a≡3(6) | a≡5(6)\n");
    printf("  --------+--------+--------+--------\n");

    for (int n_mod = 0; n_mod < 3; n_mod++) {
        printf("     %d    |", n_mod);
        for (int a_mod = 1; a_mod <= 5; a_mod += 2) {
            /* Test with actual small values */
            int n = n_mod + 3;  /* concrete n ≡ n_mod (mod 3) */
            int a = a_mod + 6;  /* concrete a ≡ a_mod (mod 6), a is odd */
            uint64_t N = 8 * n + 3;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            int c_mod3 = candidate % 3;
            printf("   %d    |", c_mod3);
        }
        printf("\n");
    }

    printf("\nInterpretation:\n");
    printf("  - When candidate mod 3 = 0, candidate is divisible by 3 (skip unless candidate = 3)\n");

    printf("\n=== Divisibility by 5 ===\n\n");
    printf("Truth table: candidate mod 5 for each (n mod 5, a mod 10):\n");
    printf("  n mod 5 | a≡1  | a≡3  | a≡5  | a≡7  | a≡9\n");
    printf("  --------+------+------+------+------+------\n");

    for (int n_mod = 0; n_mod < 5; n_mod++) {
        printf("     %d    |", n_mod);
        for (int a_mod = 1; a_mod <= 9; a_mod += 2) {
            int n = n_mod + 5;
            int a = a_mod + 10;
            uint64_t N = 8 * n + 3;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            int c_mod5 = candidate % 5;
            printf("  %d  |", c_mod5);
        }
        printf("\n");
    }

    printf("\n=== Divisibility by 7 ===\n\n");
    printf("Truth table: candidate mod 7 for each (n mod 7, a mod 14):\n");
    printf("  n mod 7 | a≡1  | a≡3  | a≡5  | a≡7  | a≡9  | a≡11 | a≡13\n");
    printf("  --------+------+------+------+------+------+------+------\n");

    for (int n_mod = 0; n_mod < 7; n_mod++) {
        printf("     %d    |", n_mod);
        for (int a_mod = 1; a_mod <= 13; a_mod += 2) {
            int n = n_mod + 7;
            int a = a_mod + 14;
            uint64_t N = 8 * n + 3;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            int c_mod7 = candidate % 7;
            printf("  %d  |", c_mod7);
        }
        printf("\n");
    }

    /* Count skippable combinations */
    printf("\n=== Skip Rate Summary ===\n\n");

    /* For mod 3: count (n,a) pairs where candidate ≡ 0 (mod 3) */
    int skip3 = 0, total3 = 0;
    for (int n_mod = 0; n_mod < 3; n_mod++) {
        for (int a_mod = 1; a_mod <= 5; a_mod += 2) {
            int n = n_mod + 3;
            int a = a_mod + 6;
            uint64_t N = 8 * n + 3;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            total3++;
            if (candidate % 3 == 0) skip3++;
        }
    }
    printf("Mod 3: %d/%d combinations have candidate ≡ 0 (mod 3) → %.1f%% skippable\n",
           skip3, total3, 100.0 * skip3 / total3);

    /* For mod 5 */
    int skip5 = 0, total5 = 0;
    for (int n_mod = 0; n_mod < 5; n_mod++) {
        for (int a_mod = 1; a_mod <= 9; a_mod += 2) {
            int n = n_mod + 5;
            int a = a_mod + 10;
            uint64_t N = 8 * n + 3;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            total5++;
            if (candidate % 5 == 0) skip5++;
        }
    }
    printf("Mod 5: %d/%d combinations have candidate ≡ 0 (mod 5) → %.1f%% skippable\n",
           skip5, total5, 100.0 * skip5 / total5);

    /* For mod 7 */
    int skip7 = 0, total7 = 0;
    for (int n_mod = 0; n_mod < 7; n_mod++) {
        for (int a_mod = 1; a_mod <= 13; a_mod += 2) {
            int n = n_mod + 7;
            int a = a_mod + 14;
            uint64_t N = 8 * n + 3;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            total7++;
            if (candidate % 7 == 0) skip7++;
        }
    }
    printf("Mod 7: %d/%d combinations have candidate ≡ 0 (mod 7) → %.1f%% skippable\n",
           skip7, total7, 100.0 * skip7 / total7);

    /* Combined skip rate for LCM(3,5,7) = 105 */
    printf("\n=== Combined Skip Rate (mod 105) ===\n\n");

    int skip_any = 0, total_combined = 0;
    for (int n_mod = 0; n_mod < 105; n_mod++) {
        for (int a_mod = 1; a_mod < 210; a_mod += 2) {  /* odd a values mod 210 */
            int n = n_mod + 105;
            int a = a_mod;
            uint64_t N = 8 * n + 3;
            if ((uint64_t)a * a >= N - 4) continue;
            uint64_t candidate = (N - (uint64_t)a * a) / 2;
            total_combined++;
            if (candidate % 3 == 0 || candidate % 5 == 0 || candidate % 7 == 0) {
                skip_any++;
            }
        }
    }
    printf("Combined: %d/%d have candidate div by 3,5, or 7 → %.1f%% skippable by wheel\n",
           skip_any, total_combined, 100.0 * skip_any / total_combined);

    printf("\nNote: This is the theoretical maximum. Actual savings depend on:\n");
    printf("  1. Cost of wheel lookup vs trial division\n");
    printf("  2. Whether we find a solution before checking all a values\n");

    return 0;
}
