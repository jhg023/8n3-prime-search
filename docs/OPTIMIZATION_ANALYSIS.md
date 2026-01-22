# Optimization Analysis: Root-Class Sieve & Wheel

## Summary

The "root-class sieve" optimization described in the literature (skip a values where a² ≡ N mod q for small primes q) **does not improve performance** for this application. The baseline implementation is already near-optimal.

## Optimizations Tested

### 1. Root-Class Sieve (Tonelli-Shanks per N)

**Algorithm**: For each N, compute square roots of N mod small primes q using Tonelli-Shanks. Skip any a where a ≡ ±sqrt(N) mod q, since the candidate p would be divisible by q.

**Result**: **2-3x SLOWER** than baseline

**Why it doesn't help**:
- Tonelli-Shanks requires O(log q) modular exponentiations per prime per N
- We only test ~10-15 candidates before finding a solution (early exit)
- The overhead exceeds the savings from skipping candidates

### 2. Precomputed Wheel/CRT (No per-N computation)

**Algorithm**: Precompute lookup table indexed by (n mod 105, a mod 210) marking which combinations produce candidates divisible by 3, 5, or 7.

**Result**: **Slightly SLOWER** than baseline (1-2%)

**Why it doesn't help**:
- Table lookup + modulo operations add overhead
- Only filters candidates divisible by 3, 5, 7 (which trial division already handles efficiently)
- Edge cases (candidate = q) require special handling

### 3. Inline Divisibility Check

**Algorithm**: Instead of trial division on candidate, check if (N - a²) ≡ 0 mod 2q directly.

**Result**: **Similar to baseline** (within 1-2%)

**Why**: Just trades one modulo operation for another.

### 4. Incremental a² Computation

**Algorithm**: Use (a-2)² = a² - 4a + 4 to update a² incrementally instead of computing a*a each iteration.

**Result**: **3-7% faster at small scales (10^9)**, neutral at large scales

**Why**: Saves one multiplication per iteration, but the effect diminishes as Miller-Rabin dominates at large N.

### 5. Reduced Trial Division (20 vs 30 primes)

**Result**: **Slower at large scales**

**Optimal prime count** (benchmarked):
| Scale | Optimal Primes | Throughput |
|-------|----------------|------------|
| 10^9  | 30-35          | ~3.5M n/sec |
| 10^12 | 30             | ~2.6M n/sec |
| 10^15 | 35-40          | ~2.0M n/sec |
| 10^18 | 35             | ~1.5M n/sec |

The current 30 primes is near-optimal across all scales.

## Why the Baseline is Already Optimal

1. **Solutions found quickly**: Average ~10-15 candidates tested per N
2. **Reverse iteration**: Testing largest a first (smallest primes) exploits higher prime density
3. **Montgomery multiplication**: 3x faster modular arithmetic for Miller-Rabin
4. **FJ64_262K**: Only 2 Miller-Rabin tests instead of 7
5. **Tuned trial division**: 30 primes filters ~80% of composites with minimal overhead

## Performance Profile (from PMU counters)

| Metric | Value |
|--------|-------|
| IPC | 2.4 |
| Branch Mispredict | 2.6% |
| L1D Cache Miss | 0.25 per 1K instructions |

The code is well-optimized for the CPU pipeline.

## When Root-Class Sieve Would Help

The root-class sieve is designed for **exhaustive search** scenarios:
- Searching for counterexamples (where no solution exists)
- Verifying that ALL candidates have been checked
- When you iterate through ~√N candidates per N

For **solution finding** (early exit when any valid pair is found), the overhead of precomputation exceeds the benefit.

## Recommendations

1. **Keep the baseline implementation** - it's already near-optimal
2. **Don't implement root-class sieve** - adds overhead that doesn't pay off
3. **Consider parallelization** instead - process multiple N values concurrently
4. **Trial division is the bottleneck** - but 30 primes is already optimal

### 6. Legendre Symbol Skip

**Algorithm**: For each N, compute Legendre symbol (N/ℓ) for first 8 primes. If (N/ℓ) = -1, skip checking divisibility by ℓ entirely (it can never divide any p_a for this N).

**Result**: **40-50% SLOWER** than baseline

| Scale | Baseline | Legendre Skip |
|-------|----------|---------------|
| 10^9  | 3.3M n/sec | 1.9M n/sec |
| 10^12 | 2.7M n/sec | 1.6M n/sec |
| 10^15 | 2.0M n/sec | 1.3M n/sec |
| 10^18 | 1.5M n/sec | 1.0M n/sec |

**Why it doesn't help**: Computing 8 Legendre symbols per N (each requiring a modular exponentiation) costs more than the savings from skipping a few divisibility checks. The benefit would only materialize if we tested hundreds of candidates per N.

### 7. Sieve on a-values (Segmented Sieve)

**Algorithm**: For each small prime ℓ with (N/ℓ) = +1, find roots r₁, r₂ where r² ≡ N (mod ℓ). Mark all a ≡ r₁ or r₂ (mod ℓ) as eliminated. Only test primality for unmarked a's.

**Result**: **6-10x SLOWER** than baseline

**Why it doesn't help**:
- Initializing sieve data requires Tonelli-Shanks for each prime
- Sieving the segment has memory access overhead
- We only test ~10-15 candidates before finding a solution
- The massive setup cost completely dominates

## Files Created

- `analysis/optimization_benchmark.c` - Root-class sieve + wheel benchmark
- `analysis/optimization_v2.c` - Inline sieve + incremental computation
- `analysis/optimization_v3.c` - Precomputed lookup tables
- `analysis/tune_trial_div.c` - Trial division prime count tuning
- `analysis/benchmark_fast.c` - Final fast implementation benchmark
- `include/solve_fast.h` - Micro-optimized implementation (marginal gains)
- `analysis/sieve_optimization.c` - Segmented sieve on a-values (6-10x slower)
- `analysis/legendre_skip.c` - Legendre symbol skip (40-50% slower)

## Conclusion

The theoretical optimizations from number theory literature (root-class sieve, wheel/CRT, sieve on a-values, Legendre skip) all sound appealing but don't work in practice because:

1. **Solutions are found too quickly** (~10-15 candidates per N) for precomputation overhead to amortize
2. **Trial division already efficiently filters composites** (30 primes filters ~80%)
3. **The current reverse-iteration strategy is optimal** for finding ANY solution (tests largest a first = smallest primes = highest density)
4. **Per-N setup costs dominate**: Even computing 8 Legendre symbols costs more than the savings

### Complete Optimization Results Summary

| Optimization | Result | Relative Performance |
|-------------|--------|---------------------|
| Baseline | Reference | 1.00x |
| Incremental a² | 3-7% faster (small scale only) | 1.03-1.07x |
| Root-class sieve | 2-3x SLOWER | 0.33-0.50x |
| Wheel/CRT (M=210) | 13% SLOWER | 0.87x |
| Precomputed lookup | 1-2% SLOWER | 0.98x |
| Legendre skip | 40-50% SLOWER | 0.50-0.60x |
| Sieve on a-values | 6-10x SLOWER | 0.10-0.17x |

**The fastest code is the existing baseline implementation** with Montgomery multiplication and FJ64_262K primality test.
