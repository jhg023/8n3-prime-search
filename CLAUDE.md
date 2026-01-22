# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

High-performance search for counterexamples to the conjecture: every integer 8n + 3 can be expressed as a² + 2p (a = odd positive integer, p = prime). Single-threaded C implementation optimized for 64-bit primality testing.

## Build Commands

```bash
make                  # Optimized release build
make debug            # Debug build with ASan/UBSan
make benchmark        # Build benchmark suite
make test             # Quick test (n = 1 to 10000)
make run-benchmark    # Run full benchmark (10M iterations/scale)
make run-benchmark-quick  # Quick benchmark (1M iterations)
make clean
```

**Running the search:**
```bash
./search                      # Default: [10^12, 10^12 + 10^7)
./search 1e9 2e9              # Custom range with scientific notation
./search 1e15 1.00001e15      # Large range
```

## Profiling (macOS)

```bash
./scripts/profile.sh [start] [end]   # Builds with -g, runs Instruments, extracts PMU data
open profile.trace                   # View in Instruments GUI
```

Uses custom Instruments templates in `~/Library/Application Support/Instruments/Templates/`:
- "My PMU Counters" (primary)
- "My CPU Counters" (fallback)

## Architecture

### Header-Only Design

All algorithmic code is in `include/` headers:
- **arith.h** - Basic 64-bit modular arithmetic (`mulmod64`, `powmod64`, `isqrt64`)
- **arith_montgomery.h** - Montgomery multiplication (3x faster modular ops for n < 2^63)
- **prime.h** - FJ64_262K primality test (2 Miller-Rabin tests via hash table witness selection)
- **solve.h** - Solution finding with multiple iteration strategies
- **fj64_table.h** - 512KB precomputed hash table for optimal MR witnesses
- **fmt.h** - Number formatting utilities

### Critical Performance Path

`find_solution()` in solve.h:
1. Iterates odd `a` values in **reverse order** (largest first = smallest prime candidates)
2. `trial_division_check()` filters ~80% of composites with 30 small primes
3. `is_prime_fj64_fast()` runs exactly 2 Montgomery-accelerated Miller-Rabin tests
4. Early exit on first valid (a, p) pair found

### Montgomery Multiplication

The `mr_witness_montgomery()` function in arith_montgomery.h:
- Uses Montgomery form for n < 2^63 (MONTGOMERY_SAFE_THRESHOLD)
- Falls back to `__uint128_t` division for n >= 2^63
- `montgomery_inverse()` computes -n^(-1) mod 2^64 via Newton's method
- `montgomery_reduce()` performs the reduction step

### FJ64_262K Algorithm

From Forisek & Jancina (2015): deterministic primality for all 64-bit integers using only 2 witnesses:
1. Always test base 2
2. Hash n to 18-bit index, lookup second witness from `fj64_bases[]` table

## Exit Codes

- `0` - No counterexamples found
- `1` - Error (invalid args, verification failure)
- `2` - Counterexample found

## Analysis Tools

The `analysis/` directory contains profiling and tuning utilities:
- `benchmark_montgomery.c` - Compare Montgomery vs standard arithmetic
- `trial_div_tuning.c` - Optimize trial division prime count
- `profile_breakdown.c` - Time breakdown by function
