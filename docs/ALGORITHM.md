# Algorithm Documentation

## Problem Statement

We search for counterexamples to the conjecture:

> For every positive integer n, there exist a positive odd integer a and a prime p such that 8n + 3 = a² + 2p.

## Search Algorithm

### Overview

For each value of n in the search range:

1. Compute N = 8n + 3
2. Iterate through odd values a in **reverse order**: √N, √N-2, ... down to 1
3. For each a, compute candidate p = (N - a²) / 2
4. Test if p ≥ 2 and p is prime
5. If any valid (a, p) pair is found, move to the next n
6. If no valid pair exists, report n as a counterexample

**Why reverse order?** Larger values of `a` produce smaller prime candidates `p`. Smaller primes are:
- Faster to verify (fewer modular operations)
- More likely to be prime (prime density ~1/ln(p))
- More likely to pass trial division quickly

### Key Observations

**Candidate form:** For odd a, we have a² ≡ 1 (mod 8). Since N ≡ 3 (mod 8):
- N - a² ≡ 2 (mod 8)
- (N - a²)/2 ≡ 1 (mod 4)

This means all prime candidates p satisfy p ≡ 1 (mod 4), which excludes p = 2 and restricts to primes of the form 4k + 1.

**Candidate range:** For the smallest a = 1:
- p = (N - 1)/2 ≈ 4n + 1 (largest candidate)

For the largest a ≈ √N:
- p ≈ 1 (smallest candidate, if valid)

### Optimizations

1. **Reverse Iteration Order**
   - Test largest `a` first (smallest prime candidates)
   - Smaller primes are faster to verify and more common
   - Typically finds solutions in fewer iterations

2. **Optimized Trial Division**
   - Test divisibility by 120 primes (3 to 661)
   - Benchmarked to find optimal prime count for large n
   - +12.4% throughput vs 30-prime baseline
   - Eliminates most composites before expensive Miller-Rabin

3. **Early Exit**
   - Stop searching for a given n once any valid (a, p) is found
   - On average, only ~10-15 candidates are tested per n

4. **FJ64_262K Primality Test**
   - Uses only 2 Miller-Rabin tests instead of 7
   - 512KB hash table selects optimal second witness
   - ~1.9x speedup over standard 7-witness test

### Trial Division Optimization

The optimal number of trial division primes was determined by benchmarking at n ~ 2×10¹⁸:

| Trial Primes | Max Prime | Throughput | vs 30 primes |
|--------------|-----------|------------|--------------|
| 30 | 127 | 504,455 n/sec | baseline |
| 50 | 233 | 530,867 n/sec | +5.2% |
| 75 | 383 | 549,413 n/sec | +8.9% |
| 100 | 547 | 557,381 n/sec | +10.5% |
| **120** | **661** | **567,220 n/sec** | **+12.4%** |
| 150 | 877 | 560,085 n/sec | +11.0% |
| 200 | 1229 | 556,434 n/sec | +10.3% |
| 308 | 2039 | 543,120 n/sec | +7.7% |

Beyond ~120 primes, the cost of additional trial divisions outweighs the benefit of filtering more composites.

## FJ64_262K Primality Test

### Background

The standard deterministic Miller-Rabin test for 64-bit integers requires 7 witnesses:
{2, 325, 9375, 28178, 450775, 9780504, 1795265022}

Forisek and Jancina (2015) showed that for any n, there exists a witness that, combined with base 2, correctly determines primality. They precomputed a hash table mapping n → optimal witness.

### Algorithm

```
FJ64_262K(n):
    if n fails Miller-Rabin test with base 2:
        return COMPOSITE
    
    h = hash(n) & 0x3FFFF        // 18-bit hash → index in [0, 262143]
    witness = table[h]           // Look up optimal witness
    
    if n fails Miller-Rabin test with witness:
        return COMPOSITE
    
    return PRIME
```

### Hash Function

```c
uint32_t hash(uint64_t x) {
    x = ((x >> 32) ^ x) * 0x45d9f3b3335b369ULL;
    x = ((x >> 32) ^ x) * 0x3335b36945d9f3bULL;
    x = ((x >> 32) ^ x);
    return x & 262143;
}
```

### Properties

- **Correctness:** Proven correct for all 64-bit integers
- **Table size:** 262,144 entries × 16 bits = 512 KB
- **Tests per candidate:** Exactly 2 (vs. up to 7 for standard test)
- **Speedup:** ~1.9x for the 8n + 3 search workload

## Miller-Rabin Test

### Single Witness Test

For odd n > 2, write n - 1 = 2^r × d where d is odd.

```
is_strong_probable_prime(n, a):
    x = a^d mod n
    
    if x == 1 or x == n - 1:
        return true
    
    for i in 1 to r - 1:
        x = x² mod n
        if x == n - 1:
            return true
        if x == 1:
            return false
    
    return false
```

### Complexity

- Modular exponentiation: O(log n) multiplications
- Each multiplication: O(1) using 128-bit arithmetic
- Total per witness: O(log n)

## Performance Analysis

### Workload Characteristics (at n ≈ 10¹²)

| Metric | Value |
|--------|-------|
| Candidates per n | ~20 |
| After trial division | ~3.5 |
| Prime ratio (of tested) | ~28.5% |

### Why FJ64_262K is Optimal

The search workload has unusually high prime density (~28.5%) because:
1. Trial division already filters most composites
2. We exit on the FIRST prime found

For workloads with many primes:
- FJ64_262K: 2 MR tests always → fast
- BPSW: 1 MR + Lucas → Lucas is expensive for primes
- MR-7: 7 MR tests → slowest

### Benchmark Results

| Algorithm | Relative Time | Tests per Prime |
|-----------|---------------|-----------------|
| FJ64_262K | 1.0x | 2 |
| BPSW | 1.2x | 1 + Lucas |
| MR-7 | 1.9x | 7 |

## References

1. Forisek, M. and Jancina, J. (2015). "Fast Primality Testing for Integers That Fit into a Machine Word." CEUR-WS Vol-1326.

2. Miller, G. L. (1976). "Riemann's hypothesis and tests for primality." Journal of Computer and System Sciences.

3. Rabin, M. O. (1980). "Probabilistic algorithm for testing primality." Journal of Number Theory.

4. Baillie, R. and Wagstaff, S. S. Jr. (1980). "Lucas Pseudoprimes." Mathematics of Computation.
