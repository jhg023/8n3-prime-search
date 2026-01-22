/*
 * SIMD/NEON Optimization Experiments
 *
 * Note: ARM NEON lacks 64-bit integer multiply (vmulq_u64),
 * so we explore alternative SIMD strategies.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <arm_neon.h>

#include "arith.h"
#include "arith_montgomery.h"
#include "prime.h"
#include "solve.h"

static inline double get_time(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

/* ========================================================================== */
/* Baseline (scalar)                                                          */
/* ========================================================================== */

static uint64_t find_solution_scalar(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    for (uint64_t a = a_max; a >= 1; a -= 2) {
        uint64_t a_sq = a * a;
        if (a_sq > N - 4) continue;
        uint64_t candidate = (N - a_sq) >> 1;
        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return a;
        }
    }
    return 0;
}

/* ========================================================================== */
/* Approach 1: Process 2 candidates with manual unrolling (no SIMD mul)       */
/* ========================================================================== */

static uint64_t find_solution_unroll2(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t limit = N - 4;
    uint64_t a = a_max;

    while (a >= 3) {
        uint64_t a0 = a, a1 = a - 2;
        uint64_t a_sq0 = a0 * a0;
        uint64_t a_sq1 = a1 * a1;

        if (a_sq0 <= limit) {
            uint64_t c0 = (N - a_sq0) >> 1;
            if (is_candidate_prime(c0)) {
                if (p_out) *p_out = c0;
                return a0;
            }
        }

        if (a_sq1 <= limit) {
            uint64_t c1 = (N - a_sq1) >> 1;
            if (is_candidate_prime(c1)) {
                if (p_out) *p_out = c1;
                return a1;
            }
        }

        a -= 4;
    }

    if (a == 1) {
        uint64_t candidate = (N - 1) >> 1;
        if (is_candidate_prime(candidate)) {
            if (p_out) *p_out = candidate;
            return 1;
        }
    }

    return 0;
}

/* ========================================================================== */
/* Approach 2: Process 4 candidates with unrolling                            */
/* ========================================================================== */

static uint64_t find_solution_unroll4(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t limit = N - 4;
    uint64_t a = a_max;

    while (a >= 7) {
        uint64_t a0 = a, a1 = a - 2, a2 = a - 4, a3 = a - 6;

        /* Compute all squares */
        uint64_t sq0 = a0 * a0, sq1 = a1 * a1, sq2 = a2 * a2, sq3 = a3 * a3;

        /* Check all 4 */
        if (sq0 <= limit) {
            uint64_t c = (N - sq0) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a0; }
        }
        if (sq1 <= limit) {
            uint64_t c = (N - sq1) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a1; }
        }
        if (sq2 <= limit) {
            uint64_t c = (N - sq2) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a2; }
        }
        if (sq3 <= limit) {
            uint64_t c = (N - sq3) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a3; }
        }

        a -= 8;
    }

    /* Handle remainder */
    while (a >= 1) {
        uint64_t a_sq = a * a;
        if (a_sq <= limit) {
            uint64_t c = (N - a_sq) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a; }
        }
        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Approach 3: SIMD for 32-bit trial division (use 4x32-bit lanes)            */
/* ========================================================================== */

/* Magic constants for divisibility (32-bit versions) */
static const uint32_t MAGIC32_3 = 0xAAAAAAABu;
static const uint32_t THRESH32_3 = 0x55555555u;
static const uint32_t MAGIC32_5 = 0xCCCCCCCDu;
static const uint32_t THRESH32_5 = 0x33333333u;
static const uint32_t MAGIC32_7 = 0x6DB6DB6Du; /* Adjusted for 32-bit */
static const uint32_t THRESH32_7 = 0x24924924u;
static const uint32_t MAGIC32_11 = 0xBA2E8BA3u;
static const uint32_t THRESH32_11 = 0x1745D174u;

/*
 * Check divisibility by 3,5,7,11 using NEON for 32-bit candidates
 * Returns bitmask: bit i set if candidate[i] passed (not divisible)
 */
static inline uint32_t simd_trial_4x32(uint32x4_t candidates) {
    /* Check div by 3 */
    uint32x4_t magic3 = vdupq_n_u32(MAGIC32_3);
    uint32x4_t thresh3 = vdupq_n_u32(THRESH32_3);
    uint32x4_t prod3 = vmulq_u32(candidates, magic3);
    uint32x4_t div3 = vcltq_u32(prod3, thresh3);

    /* Check div by 5 */
    uint32x4_t magic5 = vdupq_n_u32(MAGIC32_5);
    uint32x4_t thresh5 = vdupq_n_u32(THRESH32_5);
    uint32x4_t prod5 = vmulq_u32(candidates, magic5);
    uint32x4_t div5 = vcltq_u32(prod5, thresh5);

    /* Check div by 7 */
    uint32x4_t magic7 = vdupq_n_u32(MAGIC32_7);
    uint32x4_t thresh7 = vdupq_n_u32(THRESH32_7);
    uint32x4_t prod7 = vmulq_u32(candidates, magic7);
    uint32x4_t div7 = vcltq_u32(prod7, thresh7);

    /* Check div by 11 */
    uint32x4_t magic11 = vdupq_n_u32(MAGIC32_11);
    uint32x4_t thresh11 = vdupq_n_u32(THRESH32_11);
    uint32x4_t prod11 = vmulq_u32(candidates, magic11);
    uint32x4_t div11 = vcltq_u32(prod11, thresh11);

    /* Combine: any divisibility = skip */
    uint32x4_t any_div = vorrq_u32(vorrq_u32(div3, div5), vorrq_u32(div7, div11));

    /* Extract results: return mask where 0 = passed, ~0 = failed */
    /* We want inverse: 1 if NOT divisible */
    uint32x4_t passed = vmvnq_u32(any_div);

    /* Narrow to get a bitmask-like result */
    uint32_t p0 = vgetq_lane_u32(passed, 0) ? 1 : 0;
    uint32_t p1 = vgetq_lane_u32(passed, 1) ? 2 : 0;
    uint32_t p2 = vgetq_lane_u32(passed, 2) ? 4 : 0;
    uint32_t p3 = vgetq_lane_u32(passed, 3) ? 8 : 0;

    return p0 | p1 | p2 | p3;
}

/*
 * Batch 4 candidates, use SIMD to pre-filter by 3,5,7,11
 * Only works well when candidates fit in 32 bits
 */
static uint64_t find_solution_simd_filter(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);
    if ((a_max & 1) == 0) a_max--;

    uint64_t limit = N - 4;
    uint64_t a = a_max;

    while (a >= 7) {
        uint64_t a0 = a, a1 = a - 2, a2 = a - 4, a3 = a - 6;
        uint64_t sq0 = a0 * a0, sq1 = a1 * a1, sq2 = a2 * a2, sq3 = a3 * a3;

        /* Compute candidates */
        uint64_t c0 = (sq0 <= limit) ? (N - sq0) >> 1 : 0;
        uint64_t c1 = (sq1 <= limit) ? (N - sq1) >> 1 : 0;
        uint64_t c2 = (sq2 <= limit) ? (N - sq2) >> 1 : 0;
        uint64_t c3 = (sq3 <= limit) ? (N - sq3) >> 1 : 0;

        /* If all candidates fit in 32 bits, use SIMD filtering */
        if (c0 <= UINT32_MAX && c1 <= UINT32_MAX && c2 <= UINT32_MAX && c3 <= UINT32_MAX) {
            uint32x4_t vc = {(uint32_t)c0, (uint32_t)c1, (uint32_t)c2, (uint32_t)c3};
            uint32_t mask = simd_trial_4x32(vc);

            /* Test candidates that passed SIMD filter */
            if ((mask & 1) && c0 != 0 && is_candidate_prime(c0)) {
                if (p_out) *p_out = c0; return a0;
            }
            if ((mask & 2) && c1 != 0 && is_candidate_prime(c1)) {
                if (p_out) *p_out = c1; return a1;
            }
            if ((mask & 4) && c2 != 0 && is_candidate_prime(c2)) {
                if (p_out) *p_out = c2; return a2;
            }
            if ((mask & 8) && c3 != 0 && is_candidate_prime(c3)) {
                if (p_out) *p_out = c3; return a3;
            }
        } else {
            /* Fallback to scalar */
            if (c0 != 0 && is_candidate_prime(c0)) { if (p_out) *p_out = c0; return a0; }
            if (c1 != 0 && is_candidate_prime(c1)) { if (p_out) *p_out = c1; return a1; }
            if (c2 != 0 && is_candidate_prime(c2)) { if (p_out) *p_out = c2; return a2; }
            if (c3 != 0 && is_candidate_prime(c3)) { if (p_out) *p_out = c3; return a3; }
        }

        a -= 8;
    }

    /* Handle remainder */
    while (a >= 1) {
        uint64_t a_sq = a * a;
        if (a_sq <= limit) {
            uint64_t c = (N - a_sq) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a; }
        }
        if (a < 3) break;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Approach 4: Incremental a² with unrolling                                  */
/* ========================================================================== */

static uint64_t find_solution_incremental(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a = isqrt64(N);
    a = (a | 1);
    if (a * a > N) a -= 2;

    uint64_t limit = N - 4;
    uint64_t a_sq = a * a;

    while (a >= 1) {
        if (a_sq <= limit) {
            uint64_t candidate = (N - a_sq) >> 1;
            if (is_candidate_prime(candidate)) {
                if (p_out) *p_out = candidate;
                return a;
            }
        }

        if (a < 3) break;

        /* Incremental: (a-2)² = a² - 4a + 4 */
        a_sq = a_sq - 4 * a + 4;
        a -= 2;
    }

    return 0;
}

/* ========================================================================== */
/* Approach 5: Incremental a² with 2x unrolling                               */
/* ========================================================================== */

static uint64_t find_solution_incremental2(uint64_t n, uint64_t* p_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a = isqrt64(N);
    a = (a | 1);
    if (a * a > N) a -= 2;

    uint64_t limit = N - 4;
    uint64_t a_sq = a * a;

    while (a >= 3) {
        /* First candidate */
        if (a_sq <= limit) {
            uint64_t c = (N - a_sq) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a; }
        }

        /* Update to a-2 */
        uint64_t a_sq_next = a_sq - 4 * a + 4;
        uint64_t a_next = a - 2;

        /* Second candidate */
        if (a_sq_next <= limit) {
            uint64_t c = (N - a_sq_next) >> 1;
            if (is_candidate_prime(c)) { if (p_out) *p_out = c; return a_next; }
        }

        /* Update to a-4: (a-4)² = (a-2)² - 4(a-2) + 4 */
        a_sq = a_sq_next - 4 * a_next + 4;
        a = a_next - 2;
    }

    /* Handle a = 1 */
    if (a == 1) {
        uint64_t c = (N - 1) >> 1;
        if (is_candidate_prime(c)) { if (p_out) *p_out = c; return 1; }
    }

    return 0;
}

/* ========================================================================== */
/* Benchmark                                                                  */
/* ========================================================================== */

typedef uint64_t (*solver_fn)(uint64_t, uint64_t*);

static void benchmark(const char* name, solver_fn fn, uint64_t start_n, uint64_t count) {
    /* Verify correctness first */
    for (uint64_t i = 0; i < 1000; i++) {
        uint64_t p1, p2;
        uint64_t a1 = find_solution_scalar(start_n + i, &p1);
        uint64_t a2 = fn(start_n + i, &p2);
        uint64_t N = 8 * (start_n + i) + 3;
        if (a1 * a1 + 2 * p1 != N || a2 * a2 + 2 * p2 != N) {
            printf("  %s: VERIFICATION FAILED at n=%llu\n", name, (unsigned long long)(start_n + i));
            return;
        }
    }

    double t0 = get_time();
    volatile uint64_t sum = 0;

    for (uint64_t i = 0; i < count; i++) {
        uint64_t p;
        sum += fn(start_n + i, &p);
    }

    double elapsed = get_time() - t0;
    printf("  %-30s %10.0f n/sec  (%.3fs)\n", name, (double)count / elapsed, elapsed);
}

int main(int argc, char** argv) {
    printf("SIMD/NEON Optimization Tests\n");
    printf("============================\n\n");

    uint64_t count = 200000;
    if (argc > 1) count = strtoull(argv[1], NULL, 10);

    uint64_t scales[] = {
        1000000000ULL,
        1000000000000ULL,
        1000000000000000ULL,
        1000000000000000000ULL,
    };

    for (size_t i = 0; i < sizeof(scales) / sizeof(scales[0]); i++) {
        printf("Benchmark at n = %.2e, count = %llu\n", (double)scales[i], (unsigned long long)count);

        benchmark("Scalar (baseline)", find_solution_scalar, scales[i], count);
        benchmark("Unroll 2x", find_solution_unroll2, scales[i], count);
        benchmark("Unroll 4x", find_solution_unroll4, scales[i], count);
        benchmark("SIMD filter (32-bit)", find_solution_simd_filter, scales[i], count);
        benchmark("Incremental a²", find_solution_incremental, scales[i], count);
        benchmark("Incremental a² + unroll 2x", find_solution_incremental2, scales[i], count);

        printf("\n");
    }

    return 0;
}
