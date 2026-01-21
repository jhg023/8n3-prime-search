# Counterexample Search: 8n + 3 = a² + 2p

A high-performance search for counterexamples to the conjecture that every integer of the form 8n + 3 can be expressed as a² + 2p, where a is a positive odd integer and p is prime.

## Background

This conjecture states that for every positive integer n, the number 8n + 3 can be written as the sum of an odd square and twice a prime. For example:
- n = 1: 8(1) + 3 = 11 = 1² + 2(5)
- n = 2: 8(2) + 3 = 19 = 3² + 2(5)
- n = 3: 8(3) + 3 = 27 = 1² + 2(13)

This program systematically searches for counterexamples by testing each value of n and attempting to find a valid (a, p) pair.

## Features

- **Optimized primality testing** using FJ64_262K algorithm (Forisek-Jancina 2015)
  - Only 2 Miller-Rabin tests per candidate (vs. 7 in standard deterministic test)
  - 1.9x faster than traditional 7-witness Miller-Rabin
  - 100% deterministic for all 64-bit integers
- **Trial division pre-filtering** eliminates ~82% of composite candidates
- **Progress reporting** with throughput statistics
- **Scientific notation support** for command-line arguments

## Building

```bash
# Standard build
make

# Optimized build (recommended)
make release

# Build with debugging symbols
make debug

# Clean build artifacts
make clean
```

Or compile directly:

```bash
gcc -O3 -march=native -flto -o search src/search.c -lm -Iinclude
```

## Usage

```bash
# Default: search n ∈ [10¹², 10¹² + 10⁶)
./search

# Custom range
./search <n_start> <n_end>

# Examples
./search 1 1000000           # Search n ∈ [1, 10⁶)
./search 1e9 2e9             # Search n ∈ [10⁹, 2×10⁹)
./search 1e12 1.001e12       # Search n ∈ [10¹², 10¹² + 10⁹)
```

## Output

```
==================================================================
     Counterexample Search: 8n + 3 = a² + 2p
     (Optimized with FJ64_262K primality test)
==================================================================

Configuration:
  Range: n in [1,000,000,000,000, 1,000,001,000,000)
  Count: 1,000,000 values
  Primality test: FJ64_262K (2 Miller-Rabin tests)

Verifying known solutions...
  n=1: N=11, given (1,5), found (1,5) ... PASS
  n=2: N=19, given (3,5), found (3,5) ... PASS
  n=3: N=27, given (1,13), found (1,13) ... PASS
  n=4: N=35, given (5,5), found (1,17) ... PASS

Starting search...

n = 1,000,000,016,384 (1.6%), rate = 458,723 n/sec, max_a = 923
n = 1,000,000,049,152 (4.9%), rate = 461,205 n/sec, max_a = 1,107
...

===================================================================
RESULTS
===================================================================

Total time:           2.1 seconds
Total throughput:     461,538 n/sec
Counterexamples:      0
Maximum a observed:   1,275 (at n = 1,000,000,084,376)
```

## Performance

| Range Start | Throughput | Notes |
|-------------|------------|-------|
| 10¹² | ~460,000 n/sec | Typical performance |
| 10¹⁴ | ~380,000 n/sec | Larger candidates |
| 10¹⁶ | ~310,000 n/sec | Near 64-bit limit |

### Primality Test Comparison

The FJ64_262K algorithm provides significant speedup over alternatives:

| Algorithm | Time (relative) | Miller-Rabin Tests | Memory |
|-----------|-----------------|-------------------|--------|
| **FJ64_262K** | **1.0x** | 2 | 512 KB |
| BPSW | 1.2x | 1 + Lucas | ~340 B |
| 7-witness MR | 1.9x | 7 | None |

## Algorithm Details

### Search Strategy

For each n, the algorithm:
1. Computes N = 8n + 3
2. Iterates through odd values a = 1, 3, 5, ... up to √N
3. For each a, computes candidate p = (N - a²) / 2
4. Tests if p is prime using FJ64_262K
5. Returns the first valid (a, p) pair found, or reports a counterexample

### FJ64_262K Primality Test

The FJ64_262K algorithm (Forisek & Jancina, 2015) uses a precomputed 512KB hash table to select optimal Miller-Rabin witnesses:

1. Always test with base 2
2. Compute hash h = hash(n) to select a second witness from the table
3. Test with the selected witness

This is proven correct for all 64-bit integers and requires only 2 modular exponentiations instead of 7.

Reference: Forisek, M. and Jancina, J. (2015). "Fast Primality Testing for Integers That Fit into a Machine Word." CEUR-WS Vol-1326.

## File Structure

```
.
├── README.md
├── LICENSE
├── Makefile
├── .gitignore
├── src/
│   └── search.c           # Main search program
├── include/
│   └── fj64_table.h       # FJ64_262K hash table (512KB)
├── benchmark/
│   └── benchmark.c        # Performance benchmarks
└── docs/
    └── ALGORITHM.md       # Detailed algorithm documentation
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

## Acknowledgments

- Michal Forisek and Jakub Jancina for the FJ64_262K primality test and precomputed hash table
- The number theory community for work on prime representation conjectures
