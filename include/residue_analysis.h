/*
 * Residue Class Analysis for 8n + 3 Prime Search
 *
 * Analyzes the distribution of solutions across residue classes mod M.
 * The goal is to identify if certain residue classes are "harder" (require
 * more candidate checks) or "easier" (find solutions quickly).
 *
 * Note: A previous attempt at mod-210 wheel optimization was 1-2% SLOWER
 * due to per-N setup costs. This analysis is for research/understanding,
 * not necessarily for optimization.
 *
 * Usage:
 *   ResidueStats stats;
 *   analyze_residue_classes(1000000, 100000, &stats);
 *   print_residue_stats(&stats);
 */

#ifndef RESIDUE_ANALYSIS_H
#define RESIDUE_ANALYSIS_H

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "arith.h"
#include "prime.h"

/* ========================================================================== */
/* Configuration                                                              */
/* ========================================================================== */

/* Primary modulus for analysis: 30 = 2*3*5 */
#define RESIDUE_MODULUS_30 30

/* Extended modulus: 210 = 2*3*5*7 (for finer analysis) */
#define RESIDUE_MODULUS_210 210

/* Maximum modulus we support */
#define RESIDUE_MAX_MODULUS 210

/* ========================================================================== */
/* Data Structures                                                            */
/* ========================================================================== */

/**
 * Statistics for a single residue class r (where n % M == r)
 */
typedef struct {
    uint64_t count;             /* Number of n values in this class analyzed */
    uint64_t total_checks;      /* Total candidate checks across all n */
    uint64_t min_checks;        /* Minimum checks to find solution */
    uint64_t max_checks;        /* Maximum checks to find solution */
    uint64_t first_a_sum;       /* Sum of first a values that found solution */
    uint64_t first_a_hist[10];  /* Histogram of first successful a (1,3,5,7,9,11,...) */
} ResidueClassStats;

/**
 * Complete statistics for all residue classes mod M
 */
typedef struct {
    uint32_t modulus;                              /* The modulus M */
    ResidueClassStats classes[RESIDUE_MAX_MODULUS]; /* Stats per residue class */
    uint64_t total_n;                              /* Total n values analyzed */
    uint64_t total_checks;                         /* Total checks across all */
} ResidueStats;

/* ========================================================================== */
/* Initialization                                                             */
/* ========================================================================== */

/**
 * Initialize residue statistics structure.
 */
static inline void residue_stats_init(ResidueStats *stats, uint32_t modulus) {
    memset(stats, 0, sizeof(ResidueStats));
    stats->modulus = modulus;

    /* Initialize min_checks to max value so we can track minimum */
    for (uint32_t r = 0; r < modulus; r++) {
        stats->classes[r].min_checks = UINT64_MAX;
    }
}

/* ========================================================================== */
/* Trial Division (copied from prime.h to avoid dependency)                   */
/* ========================================================================== */

static const uint32_t RESIDUE_TRIAL_PRIMES[] = {
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
    37, 41, 43, 47, 53, 59, 61, 67, 71, 73,
    79, 83, 89, 97, 101, 103, 107, 109, 113, 127
};

static inline int residue_trial_division_check(uint64_t candidate) {
    for (int i = 0; i < 30; i++) {
        if (candidate % RESIDUE_TRIAL_PRIMES[i] == 0) {
            return (candidate == RESIDUE_TRIAL_PRIMES[i]) ? 1 : 0;
        }
    }
    return 2;
}

static inline bool residue_is_candidate_prime(uint64_t candidate) {
    int td = residue_trial_division_check(candidate);
    if (td == 0) return false;
    if (td == 1) return true;
    if (candidate <= 127) return true;
    return is_prime_fj64_fast(candidate);
}

/* ========================================================================== */
/* Analysis Functions                                                         */
/* ========================================================================== */

/**
 * Find solution for a single n and count checks.
 * Returns the number of checks made (a values tested).
 */
static inline uint64_t find_solution_count_checks(uint64_t n, uint64_t *a_out) {
    uint64_t N = 8 * n + 3;
    uint64_t a_max = isqrt64(N);

    /* Ensure a_max is odd */
    if ((a_max & 1) == 0) a_max--;

    uint64_t checks = 0;
    uint64_t a = a_max;
    uint64_t candidate = (N - a * a) >> 1;
    uint64_t delta = 2 * (a - 1);

    while (1) {
        if (candidate >= 2) {
            checks++;

            if (residue_is_candidate_prime(candidate)) {
                if (a_out) *a_out = a;
                return checks;
            }
        }

        if (a < 3) break;
        candidate += delta;
        delta -= 4;
        a -= 2;
    }

    /* Should not reach here (no counterexamples expected) */
    if (a_out) *a_out = 0;
    return checks;
}

/**
 * Analyze residue classes over a range of n values.
 * Updates stats with the results.
 */
static inline void analyze_residue_classes(uint64_t n_start, uint64_t count,
                                           ResidueStats *stats) {
    uint32_t M = stats->modulus;

    for (uint64_t i = 0; i < count; i++) {
        uint64_t n = n_start + i;
        uint32_t r = n % M;

        uint64_t a_found;
        uint64_t checks = find_solution_count_checks(n, &a_found);

        ResidueClassStats *cls = &stats->classes[r];
        cls->count++;
        cls->total_checks += checks;

        if (checks < cls->min_checks) cls->min_checks = checks;
        if (checks > cls->max_checks) cls->max_checks = checks;

        /* Track first successful a */
        if (a_found > 0) {
            cls->first_a_sum += a_found;

            /* Histogram bucket: a=1,3,5,7,... maps to bucket 0,1,2,3,... */
            uint64_t bucket = (a_found - 1) / 2;
            if (bucket < 10) {
                cls->first_a_hist[bucket]++;
            }
        }

        stats->total_n++;
        stats->total_checks += checks;
    }
}

/* ========================================================================== */
/* Reporting                                                                  */
/* ========================================================================== */

/**
 * Print summary of residue class statistics.
 */
static inline void print_residue_stats(const ResidueStats *stats) {
    uint32_t M = stats->modulus;

    printf("Residue Class Analysis (mod %u)\n", M);
    printf("=================================================\n");
    printf("Total n values analyzed: %llu\n", (unsigned long long)stats->total_n);
    printf("Total checks: %llu\n", (unsigned long long)stats->total_checks);
    printf("Global avg checks/n: %.2f\n\n",
           (double)stats->total_checks / stats->total_n);

    /* Collect active residue classes */
    typedef struct {
        uint32_t r;
        double avg_checks;
    } ClassAvg;

    ClassAvg class_avgs[RESIDUE_MAX_MODULUS];
    int num_active = 0;

    for (uint32_t r = 0; r < M; r++) {
        if (stats->classes[r].count > 0) {
            class_avgs[num_active].r = r;
            class_avgs[num_active].avg_checks =
                (double)stats->classes[r].total_checks / stats->classes[r].count;
            num_active++;
        }
    }

    /* Sort by avg_checks descending to find "hardest" classes */
    for (int i = 0; i < num_active - 1; i++) {
        for (int j = i + 1; j < num_active; j++) {
            if (class_avgs[j].avg_checks > class_avgs[i].avg_checks) {
                ClassAvg tmp = class_avgs[i];
                class_avgs[i] = class_avgs[j];
                class_avgs[j] = tmp;
            }
        }
    }

    /* Print top 10 hardest and easiest */
    printf("Top 10 HARDEST residue classes (most checks):\n");
    printf("  Class  Count       Avg Checks  Min   Max\n");
    for (int i = 0; i < 10 && i < num_active; i++) {
        uint32_t r = class_avgs[i].r;
        const ResidueClassStats *cls = &stats->classes[r];
        printf("  %3u    %-10llu  %-10.2f  %-5llu %llu\n",
               r, (unsigned long long)cls->count, class_avgs[i].avg_checks,
               (unsigned long long)cls->min_checks,
               (unsigned long long)cls->max_checks);
    }

    printf("\nTop 10 EASIEST residue classes (fewest checks):\n");
    printf("  Class  Count       Avg Checks  Min   Max\n");
    for (int i = num_active - 1; i >= 0 && i >= num_active - 10; i--) {
        uint32_t r = class_avgs[i].r;
        const ResidueClassStats *cls = &stats->classes[r];
        printf("  %3u    %-10llu  %-10.2f  %-5llu %llu\n",
               r, (unsigned long long)cls->count, class_avgs[i].avg_checks,
               (unsigned long long)cls->min_checks,
               (unsigned long long)cls->max_checks);
    }

    /* Print overall statistics */
    double global_avg = (double)stats->total_checks / stats->total_n;
    double hardest_avg = class_avgs[0].avg_checks;
    double easiest_avg = class_avgs[num_active - 1].avg_checks;

    printf("\nSummary:\n");
    printf("  Active residue classes: %d / %u\n", num_active, M);
    printf("  Hardest class avg: %.2f checks (%.1f%% above global)\n",
           hardest_avg, 100.0 * (hardest_avg - global_avg) / global_avg);
    printf("  Easiest class avg: %.2f checks (%.1f%% below global)\n",
           easiest_avg, 100.0 * (global_avg - easiest_avg) / global_avg);
}

/**
 * Print detailed histogram of first successful a values.
 */
static inline void print_first_a_histogram(const ResidueStats *stats) {
    printf("\nFirst Successful 'a' Value Distribution:\n");
    printf("=========================================\n");

    /* Aggregate histogram across all classes */
    uint64_t total_hist[10] = {0};
    uint64_t total_count = 0;

    for (uint32_t r = 0; r < stats->modulus; r++) {
        for (int i = 0; i < 10; i++) {
            total_hist[i] += stats->classes[r].first_a_hist[i];
            total_count += stats->classes[r].first_a_hist[i];
        }
    }

    if (total_count == 0) {
        printf("  (no data)\n");
        return;
    }

    printf("  a value  Count        Percentage\n");
    for (int i = 0; i < 10; i++) {
        uint64_t a = 2 * i + 1;  /* a = 1, 3, 5, 7, 9, 11, 13, 15, 17, 19 */
        double pct = 100.0 * total_hist[i] / total_count;
        printf("  a=%-3llu    %-12llu %.2f%%\n",
               (unsigned long long)a, (unsigned long long)total_hist[i], pct);
    }

    printf("\nNote: Largest 'a' values (first checked) produce smallest 'p' candidates.\n");
    printf("High counts for large 'a' mean small primes are commonly solutions.\n");
}

#endif /* RESIDUE_ANALYSIS_H */
