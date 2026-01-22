/*
 * Compare trial division loop unrolling
 */

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

#include "../include/arith.h"
#include "../include/arith_montgomery.h"
#include "../include/fj64_table.h"

/* Trial primes */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};
#define NUM_TRIAL_PRIMES 30

#define ITERATIONS 10000000ULL

static inline uint32_t fj64_hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}

static inline bool is_prime_fj64(uint64_t n) {
    uint64_t n_inv = montgomery_inverse(n);
    uint64_t r_sq = montgomery_r_squared(n);
    if (!mr_witness_montgomery_cached(n, 2, n_inv, r_sq))
        return false;
    return mr_witness_montgomery_cached(n, fj64_bases[fj64_hash(n)], n_inv, r_sq);
}

/* Version A: Current (inline 3, loop rest) */
static inline int trial_division_current(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;
    for (int i = 3; i < NUM_TRIAL_PRIMES; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* Version B: Unrolled 4x */
static inline int trial_division_unrolled(uint64_t candidate) {
    if (candidate % 3 == 0) return (candidate == 3) ? 1 : 0;
    if (candidate % 5 == 0) return (candidate == 5) ? 1 : 0;
    if (candidate % 7 == 0) return (candidate == 7) ? 1 : 0;

    int i = 3;
    for (; i < 27; i += 4) {
        if (candidate % TRIAL_PRIMES[i] == 0)
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+1] == 0)
            return (candidate == TRIAL_PRIMES[i+1]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+2] == 0)
            return (candidate == TRIAL_PRIMES[i+2]) ? 1 : 0;
        if (candidate % TRIAL_PRIMES[i+3] == 0)
            return (candidate == TRIAL_PRIMES[i+3]) ? 1 : 0;
    }
    if (candidate % TRIAL_PRIMES[27] == 0)
        return (candidate == TRIAL_PRIMES[27]) ? 1 : 0;
    if (candidate % TRIAL_PRIMES[28] == 0)
        return (candidate == TRIAL_PRIMES[28]) ? 1 : 0;
    if (candidate % TRIAL_PRIMES[29] == 0)
        return (candidate == TRIAL_PRIMES[29]) ? 1 : 0;
    return 2;
}

static inline bool is_candidate_prime_current(uint64_t candidate) {
    int td = trial_division_current(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64(candidate);
}

static inline bool is_candidate_prime_unrolled(uint64_t candidate) {
    int td = trial_division_unrolled(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64(candidate);
}

static inline uint64_t find_solution_current(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_current(candidate)) {
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

static inline uint64_t find_solution_unrolled(uint64_t N, uint64_t a_max, uint64_t* p_out) {
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            if (is_candidate_prime_unrolled(candidate)) {
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

double benchmark(uint64_t n_start, uint64_t (*find_fn)(uint64_t, uint64_t, uint64_t*)) {
    clock_t start = clock();
    volatile uint64_t sum = 0;

    uint64_t N = 8 * n_start + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t i = 0; i < ITERATIONS; i++) {
        uint64_t p;
        sum += find_fn(N, a_max, &p);

        N += 8;
        uint64_t next_a = a_max + 2;
        if (next_a * next_a <= N) a_max = next_a;
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    return ITERATIONS / elapsed;
}

int main(void) {
    printf("Comparing trial division loop unrolling (%llu iterations)\n\n", ITERATIONS);

    uint64_t scales[] = {1000000000000ULL, 1000000000000000ULL, 2000000000000000000ULL};
    const char* names[] = {"10^12", "10^15", "2e18"};

    for (int s = 0; s < 3; s++) {
        printf("=== n = %s ===\n", names[s]);

        double current_sum = 0, unrolled_sum = 0;
        int runs = 5;

        for (int i = 0; i < runs; i++) {
            double c = benchmark(scales[s], find_solution_current);
            double u = benchmark(scales[s], find_solution_unrolled);
            current_sum += c;
            unrolled_sum += u;
            printf("  Run %d: current=%.0f, unrolled=%.0f (%.2f%%)\n",
                   i+1, c, u, 100.0*(u-c)/c);
        }

        double current_avg = current_sum / runs;
        double unrolled_avg = unrolled_sum / runs;
        printf("  Average: current=%.0f, unrolled=%.0f (%.2f%%)\n\n",
               current_avg, unrolled_avg, 100.0*(unrolled_avg-current_avg)/current_avg);
    }

    return 0;
}
