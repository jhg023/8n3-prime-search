# Counterexample Search: 8n + 3 = a^2 + 2p

A high-performance search for counterexamples to the conjecture that every integer of the form 8n + 3 can be expressed as a^2 + 2p, where a is a positive odd integer and p is prime.

## Background

This conjecture states that for every positive integer n, the number 8n + 3 can be written as the sum of an odd square and twice a prime. For example:
- n = 1: 8(1) + 3 = 11 = 1^2 + 2(5)
- n = 2: 8(2) + 3 = 19 = 3^2 + 2(5)
- n = 3: 8(3) + 3 = 27 = 1^2 + 2(13)

This program systematically searches for counterexamples by testing each value of n and attempting to find a valid (a, p) pair.

## Features

- **Optimized primality testing** using FJ64_262K algorithm (Forisek-Jancina 2015)
  - Only 2 Miller-Rabin tests per candidate (vs. 7 in standard deterministic test)
  - 100% deterministic for all 64-bit integers
- **Montgomery multiplication** for 3x faster modular arithmetic (n < 2^63)
  - Hybrid implementation falls back to `__uint128_t` for larger moduli
- **Optimized trial division** with 120 primes (up to 661)
- **Reverse iteration** tests smallest prime candidates first for faster solutions
- **Progress reporting** with throughput and 32-bit candidate statistics
- **Scientific notation support** for command-line arguments
- **Benchmark suite** for comparing performance across scales

## Building

```bash
# Build optimized release (default)
make

# Build with debugging symbols and sanitizers
make debug

# Build benchmark suite
make benchmark

# Clean build artifacts
make clean
```

## Usage

```bash
# Default: search n in [10^12, 10^12 + 10^7)
./search

# Custom range
./search <n_start> <n_end>

# Examples
./search 1 1000000           # Search n in [1, 10^6)
./search 1e9 2e9             # Search n in [10^9, 2*10^9)
./search 1e15 1.00001e15     # Search n in [10^15, 10^15 + 10^10)
```

## Benchmarking

The benchmark suite tests throughput at various scales from 10^6 to 2*10^18:

```bash
# Run benchmark (10M iterations per scale, ~1 minute)
make run-benchmark

# Quick benchmark (1M iterations per scale, ~10 seconds)
make run-benchmark-quick

# Custom iteration count
./benchmark/benchmark_suite --count 5000000
```

Sample output:
```
Benchmark: 8n + 3 = a^2 + 2p
Iterations per scale: 1,000,000

Scale       Bits     Rate (n/sec)    Avg checks  Time (s)
--------------------------------------------------------------
10^6          23        4,474,433          5.86      0.22
10^9          33        3,211,180          8.22      0.31
10^12         43        2,285,129         10.36      0.44
10^15         53        1,788,249         12.97      0.56
10^17         60        1,596,141         13.85      0.63
10^18         63        1,388,892         15.51      0.72
2e18          64        1,361,168         15.65      0.73
--------------------------------------------------------------
```

## Performance

| Range Start | Throughput | Notes |
|-------------|------------|-------|
| 10^6 | ~4,500,000 n/sec | Small numbers, fast |
| 10^12 | ~2,300,000 n/sec | All candidates fit in 32-bit |
| 10^15 | ~1,800,000 n/sec | Mixed 32/64-bit candidates |
| 2*10^18 | ~1,360,000 n/sec | Near 64-bit limit |

## Algorithm Details

### Search Strategy

For each n, the algorithm:
1. Computes N = 8n + 3
2. Iterates through odd values a in **reverse order** (from sqrt(N) down to 1)
3. For each a, computes candidate p = (N - a^2) / 2
4. Applies trial division with 120 small primes to filter composites
5. Tests remaining candidates with FJ64_262K primality test
6. Returns the first valid (a, p) pair found, or reports a counterexample

The reverse iteration order tests smaller prime candidates first, which are faster to verify and more likely to be prime.

### FJ64_262K Primality Test

The FJ64_262K algorithm (Forisek & Jancina, 2015) uses a precomputed 512KB hash table to select optimal Miller-Rabin witnesses:

1. Always test with base 2
2. Compute hash h = hash(n) to select a second witness from the table
3. Test with the selected witness

This is proven correct for all 64-bit integers and requires only 2 modular exponentiations instead of 7.

Reference: Forisek, M. and Jancina, J. (2015). "Fast Primality Testing for Integers That Fit into a Machine Word." CEUR-WS Vol-1326.

### Montgomery Multiplication

For moduli n < 2^63, the Miller-Rabin tests use Montgomery multiplication instead of standard `__uint128_t` division. This replaces expensive division operations (~17 cycles) with multiplication and shifts (~9 cycles), providing approximately 3x faster modular arithmetic.

For moduli n ≥ 2^63, the implementation falls back to standard `__uint128_t` arithmetic to avoid overflow in the Montgomery reduction step.

## File Structure

```
.
├── README.md
├── LICENSE
├── Makefile
├── CHANGELOG.md
├── src/
│   └── search.c              # Main search program
├── include/
│   ├── arith.h               # Arithmetic utilities (mulmod, powmod, isqrt)
│   ├── arith_montgomery.h    # Montgomery multiplication (3x faster modular ops)
│   ├── prime.h               # Primality testing (FJ64_262K, trial division)
│   ├── solve.h               # Solution finding strategies
│   ├── fmt.h                 # Number formatting utilities
│   └── fj64_table.h          # FJ64_262K hash table (512KB)
├── benchmark/
│   └── benchmark_suite.c     # Performance benchmarks
├── analysis/
│   ├── prime_sizes.c         # Prime candidate size analysis
│   ├── profile_breakdown.c   # Time breakdown profiler
│   ├── trial_div_profile.c   # Trial division hit rate analysis
│   ├── wheel_analysis.c      # Wheel factorization potential analysis
│   └── benchmark_montgomery.c # Montgomery vs standard comparison
└── docs/
    └── ALGORITHM.md          # Detailed algorithm documentation
```

## Exit Codes

- `0` - Search completed, no counterexamples found
- `1` - Error (invalid arguments, verification failure)
- `2` - Counterexample found!

## Contributing

Contributions are welcome! Areas of interest:
- Multi-threaded parallelization
- GPU acceleration (CUDA/OpenCL)
- Extended range searching (beyond 64-bit)
- Alternative primality tests

## License

MIT License - see [LICENSE](LICENSE) for details.

## References

1. Forisek, M. and Jancina, J. (2015). "Fast Primality Testing for Integers That Fit into a Machine Word." CEUR-WS Vol-1326. http://ceur-ws.org/Vol-1326/020-Forisek.pdf

2. Original conjecture and related work on representing integers as sums of squares and primes.
