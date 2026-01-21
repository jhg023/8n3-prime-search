# Optimization Analysis for 64-bit Search

## Current Performance Profile (at n ~ 2e18)

| Component | Time % | Notes |
|-----------|--------|-------|
| Miller-Rabin | 71.2% | 2 tests per candidate that passes trial div |
| Trial Division | 16.4% | 120 primes, filters 86.5% of candidates |
| Iteration | 12.4% | Loop overhead, isqrt, candidate computation |

**Throughput**: ~590k n/sec at 64-bit, ~15.75 checks per n

## Optimization Options

### 1. Wheel Factorization (Skip certain a values)

**What it does**: Pre-compute which (n mod 105, a mod 210) pairs produce candidates
divisible by 3, 5, or 7. Skip those a values entirely.

**Theoretical skip rate**: 54.5% of a-values

**Estimated speedup**:
- Saves: 54.5% of (12.4% + 16.4%) = 15.7% of total time
- New throughput: ~700k n/sec (+18%)

**Complexity**: Medium (precompute lookup table, modify iteration logic)

**Caveat**: We find solutions quickly (avg 15.75 checks), so we rarely iterate
through all a values. The wheel helps most when searching for counterexamples
(which don't exist).

### 2. Reduce Trial Division Primes

**Current**: 120 primes filter 86.5% of composites

**Analysis**: The marginal benefit of later primes is low. Primes 3,5,7,11,13
alone filter ~77% of composites. Primes beyond 127 rarely catch anything new.

**Experiment needed**: Profile hit rate of each trial prime range.

**Estimated speedup**: 5-10% if we can reduce to ~50 primes

### 3. Montgomery Multiplication

**Current**: `mulmod64` uses `__uint128_t` division (~15-20 cycles per op)

**Montgomery**: Avoids division, uses multiplication + shifts (~8-10 cycles)

**Impact on Miller-Rabin**: Each MR test does ~60 modular multiplications.
  - Current: ~60 × 17 = 1020 cycles per MR
  - Montgomery: ~60 × 9 + setup cost = 540 + 100 = 640 cycles per MR

**Estimated speedup on MR**: ~37%

**Overall speedup**: 37% × 71.2% = 26% faster

**Complexity**: High (requires implementing Montgomery reduction correctly)

### 4. Batch Sieving for Small Primes

**Idea**: For consecutive n values, the candidates form arithmetic progressions.
We could sieve a range of candidates to identify those divisible by small primes
without testing each one individually.

**Challenge**: Access pattern is scattered. Each n produces multiple candidates
at different locations.

**Estimated speedup**: Unknown, likely complex to implement correctly

### 5. Better Early Termination

**Current**: We iterate from largest a (smallest candidate) until we find a prime.

**Observation**: At n ~ 2e18, we check ~15.75 candidates before finding a solution.
This is already good - we're leveraging the higher prime density at small numbers.

**Potential**: Not much room for improvement here.

## Recommended Priority

1. **Montgomery multiplication** (26% estimated speedup, high effort)
   - Biggest impact because MR dominates
   - Well-understood algorithm

2. **Wheel factorization** (18% estimated speedup, medium effort)
   - Reduces wasted work
   - Simple lookup table

3. **Trial division tuning** (5-10% speedup, low effort)
   - Profile and optimize prime count
   - Low risk

## Combined Estimate

If we implement Montgomery + Wheel:
- Montgomery: 26% speedup on MR → total 18.5% faster
- Wheel: 15.7% of remaining → total 13% faster

Combined: ~35% speedup → ~800k n/sec at 64-bit (from 590k)

## Quick Win: Trial Division Profiling

Let's first verify whether our 120 primes are optimal or if we're doing
unnecessary work.
