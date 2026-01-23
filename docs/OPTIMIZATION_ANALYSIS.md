# Optimization Analysis: 8n+3 Prime Search

This document records all optimization experiments, whether successful or not. The purpose is to prevent re-trying approaches that have already been tested.

## Summary

The baseline implementation with Montgomery multiplication and FJ64_262K primality testing is highly optimized. Most theoretical optimizations from number theory literature **do not help** because solutions are found too quickly (~10-15 candidates per N) for precomputation to amortize.

## Current Performance (v2.0)

| Scale | Throughput | Avg Candidates/N |
|-------|------------|------------------|
| 10^6  | ~5.8M n/sec | ~6 |
| 10^12 | ~2.8M n/sec | ~12 |
| 10^15 | ~2.1M n/sec | ~15 |
| 10^18 | ~1.6M n/sec | ~16 |

---

## FAILED OPTIMIZATIONS (Do Not Retry)

### 1. Root-Class Sieve (Tonelli-Shanks per N)

**Algorithm**: For each N, compute square roots of N mod small primes q using Tonelli-Shanks. Skip any a where a ≡ ±sqrt(N) mod q, since the candidate p would be divisible by q.

**Result**: **2-3x SLOWER** than baseline

**Why it fails**:
- Tonelli-Shanks requires O(log q) modular exponentiations per prime per N
- We only test ~10-15 candidates before finding a solution (early exit)
- Setup cost completely dominates any savings

### 2. Precomputed Wheel/CRT (M=210)

**Algorithm**: Precompute lookup table indexed by (n mod 105, a mod 210) marking which combinations produce candidates divisible by 3, 5, or 7.

**Result**: **1-2% SLOWER** than baseline (13% slower in some tests)

**Why it fails**:
- Table lookup + modulo operations add overhead
- Only filters candidates divisible by 3, 5, 7 (which trial division already handles efficiently)
- Edge cases (candidate = q) require special handling

### 3. Legendre Symbol Skip

**Algorithm**: For each N, compute Legendre symbol (N/ℓ) for first 8 primes. If (N/ℓ) = -1, skip checking divisibility by ℓ.

**Result**: **40-50% SLOWER** than baseline

| Scale | Baseline | Legendre Skip |
|-------|----------|---------------|
| 10^9  | 3.3M n/sec | 1.9M n/sec |
| 10^12 | 2.7M n/sec | 1.6M n/sec |
| 10^15 | 2.0M n/sec | 1.3M n/sec |
| 10^18 | 1.5M n/sec | 1.0M n/sec |

**Why it fails**: Computing 8 Legendre symbols per N (each requiring a modular exponentiation) costs more than the savings from skipping divisibility checks. Benefit only materializes if testing hundreds of candidates per N.

### 4. Sieve on a-values (Segmented Sieve)

**Algorithm**: For each small prime ℓ with (N/ℓ) = +1, find roots r₁, r₂ where r² ≡ N (mod ℓ) using Tonelli-Shanks. Mark all a ≡ r₁ or r₂ (mod ℓ) as eliminated. Only test primality for unmarked a's.

**Result**: **6-10x SLOWER** than baseline

**Why it fails**:
- Initializing sieve data requires Tonelli-Shanks for each prime
- Sieving the segment has memory access overhead
- We only test ~10-15 candidates before finding a solution
- Massive setup cost completely dominates

### 5. Bitset Sieve (Batch K Candidates)

**Algorithm**: For K=32-64 candidates, precompute a bitmask of survivors (candidates not divisible by any small prime) using parallel trial division, then only run primality tests on survivors.

**Result**: **No improvement** or slightly slower

**Why it fails**:
- Overhead of building survivor mask exceeds benefit
- We find solutions quickly, so rarely process full batches
- Sequential early-exit trial division is already efficient

### 6. Incremental Residue Updates for Trial Division

**Algorithm**: Instead of computing candidate % prime, maintain residues incrementally as a decreases by 2. Uses addition/subtraction instead of modulo.

**Result**: **36% SLOWER** than baseline

**Why it fails**:
- Updating 30 residues per iteration (additions + conditional subtractions) costs more than the modulo operations it replaces
- Modern CPUs have efficient division hardware
- Compiler already optimizes small constant divisors to multiply+shift

### 7. libdivide for Trial Division

**Algorithm**: Use libdivide library to replace modulo with multiply+shift operations.

**Result**: **No significant improvement**

**Why it fails**:
- Compiler already applies similar optimizations for constant divisors
- For 64-bit dividends with small constant divisors, hardware is competitive
- libdivide overhead (struct lookup) offsets any gains

### 8. GCD with Primorial

**Algorithm**: Check `gcd(candidate, 223092870)` (product of first 8 primes) first. If GCD > 1, determine which prime divides it.

**Result**: **Slower** than direct trial division

**Why it fails**:
- GCD algorithm requires multiple iterations
- Still need to identify which prime after GCD
- Direct modulo checks with early exit are faster

### 9. Extended Trial Division (120 primes up to 661)

**Algorithm**: Use 120 primes instead of 30 for trial division.

**Result**: **Slower at large scales**

| Primes | 10^9 | 10^12 | 10^15 | 10^18 |
|--------|------|-------|-------|-------|
| 30     | Best | Best  | Best  | Best  |
| 120    | OK   | Slower| Slower| Slower|

**Why it fails**: Additional primes rarely catch composites that pass the first 30. The marginal benefit decreases rapidly while overhead increases linearly.

### 10. Only 8 Trial Primes

**Algorithm**: Use only first 8 primes (3-23) for trial division, filtering ~75% of composites.

**Result**: **Slower at large scales**

**Why it fails**: Too many composites pass to Miller-Rabin, which is expensive. 30 primes is the optimal balance.

### 11. SIMD/NEON Optimization

**Algorithm**: Use ARM NEON vector instructions to process multiple candidates in parallel.

**Result**: **Not beneficial**

**Why it fails**:
- ARM NEON lacks 64-bit integer multiply (`vmulq_u64`)
- Miller-Rabin is inherently sequential (modular exponentiation)
- Manual unrolling (2x, 4x) shows no improvement over scalar code
- Compiler already vectorizes where possible

### 12. Simplified FJ64 Hash

**Algorithm**: Reduce FJ64 hash from 3 rounds to 1-2 rounds.

**Result**: **Cannot use** - Would break correctness

The FJ64_262K hash is carefully designed to work with the precomputed witness table. Changing the hash function would require regenerating the table.

---

## SUCCESSFUL OPTIMIZATIONS (Already Implemented)

### 1. Montgomery Multiplication (v1.0)
- Replaces `__uint128_t` division with multiply+shift
- **~3x faster** modular arithmetic for Miller-Rabin
- 37% reduction in MR time, ~26% overall speedup

### 2. FJ64_262K Primality Test (v1.0)
- Only 2 Miller-Rabin tests instead of 7
- Hash-based witness selection from 512KB table
- Deterministic for all 64-bit integers

### 3. Trial Division with 30 Primes (v1.0)
- Optimal balance tested across scales
- Filters ~80% of composites
- Primes 3-127

### 4. Reverse Iteration (Largest a First) (v1.0)
- Tests largest a values first (smallest candidate primes)
- Exploits higher prime density at small numbers
- Finds solutions in ~10-15 candidates on average

### 5. Montgomery Constant Caching (v1.3.1)
- Cache `n_inv` and `r_sq` across both MR witnesses
- Saves ~8-12 ns per candidate

### 6. Incremental N and a_max Tracking (v1.3.2)
- Maintain N = 8n + 3 and a_max incrementally
- Avoid recomputing isqrt64(N) each iteration
- **3-7% speedup** depending on scale

### 7. Branchless MR Exponentiation (v1.3.3)
- Replace conditional branch with conditional select
- Avoids branch mispredictions on unpredictable exponent bits
- **~3-4% speedup** at large scales

### 8. Inline First 3-7 Trial Primes (v1.3.4)
- Inline divisibility checks for 3, 5, 7 (catch ~50% of composites)
- Avoids loop overhead for most common filtering
- **~2% speedup**

### 9. Trial Division Loop Unrolling (v1.3.4)
- Unroll remaining trial division 4x
- Reduces branch overhead
- **~1% speedup**

### 10. FJ64 Hash Prefetch (v1.3.4)
- Prefetch FJ64 hash table entry while computing Montgomery constants
- Hides memory latency by overlapping cache fetch with CPU work
- **~1% speedup**

### 11. OpenMP Parallelization (v2.0)
- Parallelize search with OpenMP for multi-core systems
- Near-linear scaling with core count
- **~10x speedup** on 14-core system

---

## MARGINAL / NEUTRAL OPTIMIZATIONS

### 1. Inline Divisibility Check
Instead of trial division on candidate, check if (N - a²) ≡ 0 mod 2q directly.

**Result**: Within 1-2% of baseline - trades one modulo for another.

### 2. Incremental a² Computation
Use (a-2)² = a² - 4a + 4 to update a² incrementally.

**Result**: **3-7% faster at small scales (10^6-10^9)**, neutral at large scales. Effect diminishes as Miller-Rabin dominates.

### 3. Skip base-2 MR when hash witness == 2
When FJ64 hash gives witness 2, only do one MR test.

**Result**: Marginal improvement. Witness == 2 occurs only ~0.4% of the time.

### 4. Deferred Hash Computation
Compute FJ64 hash only after trial division passes.

**Result**: Neutral. Hash computation is cheap and already overlapped with prefetch.

### 5. Fast 32-bit Montgomery Inverse
Use 4 Newton iterations instead of 5 for candidates < 2^32.

**Result**: Neutral. Most candidates at large scales are > 2^32.

### 6. Progress Check Interval Tuning
Tested intervals: 16K, 64K, 256K, 1M iterations.

**Result**: 256K (0x3FFFF) chosen. Overhead is <0.5% regardless of interval.

### 7. Removing candidate >= 2 Check
At n >= 10^9, minimum candidate is always > 2.

**Result**: Negligible improvement. Check is cheap and branch predictor handles it.

---

## PERFORMANCE PROFILE

From PMU counters at n ~ 10^12:

| Metric | Value |
|--------|-------|
| IPC | 2.4 |
| Branch Mispredict | 2.6% |
| L1D Cache Miss | 0.25 per 1K instructions |

Time breakdown:
| Component | Percentage |
|-----------|------------|
| Miller-Rabin | ~50% |
| Trial Division | ~25% |
| Iteration Overhead | ~25% |

---

## WHY THEORETICAL OPTIMIZATIONS FAIL

The root-class sieve and related number-theoretic optimizations are designed for **exhaustive search** scenarios:
- Searching for counterexamples (where no solution exists)
- Verifying that ALL candidates have been checked
- When you iterate through ~√N candidates per N

For **solution finding** with early exit (our use case):
1. Solutions are found too quickly (~10-15 candidates per N)
2. Precomputation cost cannot amortize
3. Setup overhead exceeds filtering benefit
4. The current reverse-iteration strategy is near-optimal

---

## RECOMMENDATIONS

1. **Keep the current implementation** - it's already near-optimal for solution finding
2. **Don't re-implement**: root-class sieve, Legendre skip, sieve on a-values, bitset sieve, incremental residues
3. **Parallelization is the only remaining lever** - process multiple N values concurrently (done in v2.0)
4. **For counterexample search** (if needed): would require different strategy, possibly the sieve approaches

---

## FILES

### Removed (findings documented above)

**Legendre/Sieve experiments (all slower):**
- `legendre_test.c`, `legendre_test2.c`, `legendre_test3.c`
- `legendre_skip.c`, `sieve_optimization.c`
- `bitset_sieve.c`, `residue_trial_div.c`

**Optimization experiments:**
- `optimization_benchmark.c`, `optimization_v2.c`, `optimization_v3.c`
- `optimization_ideas.c`, `optimization_ideas2.c`, `optimization_ideas3.c`
- `optimization_test.c`, `advanced_opt_test.c`, `final_optimization.c`
- `compare_deferred_hash.c`, `compare_unroll.c`

**Trial division tuning (conclusion: 30 primes optimal):**
- `trial_div_final.c`, `trial_div_focused.c`, `trial_div_profile.c`
- `trial_div_hitrate.c`, `tune_trial_div.c`

**Other one-off tests:**
- `progress_overhead.c` through `progress_overhead4.c`
- `libdivide_trial_div.c`, `libdivide_full_search.c`
- `simd_test.c`, `asm_check.c`, `line_profile.c`
- `find_solution_overhead.c`, `isqrt_overhead.c`
- `prime_sizes.c`, `wheel_analysis.c`

**Unused code:**
- `include/solve_fast.h`, `analysis/benchmark_fast.c`
- `analysis/optimization_estimate.md`

### Retained in analysis/

- `benchmark_montgomery.c` - Montgomery vs standard arithmetic comparison
- `profile_breakdown.c` - Time breakdown by component
- `trial_div_tuning.c` - Representative trial division tuning (30 primes conclusion)
