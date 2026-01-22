/*
 * Minimal file to inspect assembly of hot functions
 */
#include <stdint.h>
#include <stdbool.h>

/* Montgomery reduce - HOT PATH */
__attribute__((noinline))
uint64_t montgomery_reduce(__uint128_t t, uint64_t n, uint64_t n_inv) {
    uint64_t m = (uint64_t)t * n_inv;
    __uint128_t mn = (__uint128_t)m * n;
    uint64_t u = (t + mn) >> 64;
    return (u >= n) ? (u - n) : u;
}

/* Montgomery multiply - HOT PATH */
__attribute__((noinline))
uint64_t montgomery_mul(uint64_t a, uint64_t b, uint64_t n, uint64_t n_inv) {
    __uint128_t t = (__uint128_t)a * b;
    return montgomery_reduce(t, n, n_inv);
}

/* Montgomery inverse */
__attribute__((noinline))
uint64_t montgomery_inverse(uint64_t n) {
    uint64_t x = n;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    x *= 2 - n * x;
    return -x;
}

/* Trial division - HOT PATH */
static const uint32_t TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};

__attribute__((noinline))
int trial_division_check(uint64_t candidate) {
    for (int i = 0; i < 30; i++) {
        if (candidate % TRIAL_PRIMES[i] == 0) {
            return (candidate == TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

/* isqrt64 */
__attribute__((noinline))
uint64_t isqrt64(uint64_t n) {
    if (n == 0) return 0;
    uint64_t x = n;
    uint64_t y = (x + 1) >> 1;
    while (y < x) {
        x = y;
        y = (x + n / x) >> 1;
    }
    return x;
}

/* Main loop structure */
__attribute__((noinline))
uint64_t find_solution_check(uint64_t n) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t a = a_max;
    while (1) {
        uint64_t a_sq = a * a;
        if (a_sq <= N - 4) {
            uint64_t candidate = (N - a_sq) >> 1;
            int td = trial_division_check(candidate);
            if (td == 1) return a;  /* Small prime */
            if (td == 2) return a;  /* Needs MR - pretend it's prime for this test */
        }
        if (a < 3) break;
        a -= 2;
    }
    return 0;
}
