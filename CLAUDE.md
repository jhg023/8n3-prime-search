# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

High-performance search for counterexamples to the conjecture: every integer 8n + 3 can be expressed as a² + 2p (a = odd positive integer, p = prime). Multi-threaded CPU (OpenMP) and GPU (Metal) implementations optimized for 64-bit primality testing.

## Build Commands

```bash
make                  # Optimized release build (CPU with OpenMP)
make metal            # GPU-accelerated build (macOS Metal)
make debug            # Debug build with ASan/UBSan
make benchmark        # Build benchmark suite
make test             # Quick test (n = 1 to 10000)
make test-gpu         # Test GPU against CPU
make run-benchmark    # Run full benchmark (10M iterations/scale)
make run-benchmark-quick  # Quick benchmark (1M iterations)
make clean
```

**Running the CPU search:**
```bash
./search                      # Default: [10^12, 10^12 + 10^7)
./search 1e9 2e9              # Custom range with scientific notation
./search 1e9 2e9 --threads 4  # Specify thread count
```

**Running the GPU search (macOS only):**
```bash
./search_gpu                        # Default range on GPU
./search_gpu 1e9 2e9                # Custom range
./search_gpu 1e9 2e9 --batch-size 100000  # Adjust batch size
./search_gpu --verify-only          # Run verification only
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

## Metal GPU Implementation

The `metal/` directory contains the GPU compute implementation:

### Files
- **prime_search.metal** - Metal compute shader with 128-bit software arithmetic
- **metal_host.m** - Objective-C Metal API interface
- **metal_host.h** - C-compatible header for GPU functions

### Architecture
- One threadgroup per `n` value, 256 threads per group
- Each thread tests different `a` values (largest-first strategy preserved)
- Atomic early termination when solution found
- All potential counterexamples verified on CPU before reporting

### Key Implementation Details
- **128-bit arithmetic**: Software-emulated using 4 32x32 multiplies (Metal lacks native 128-bit)
- **Montgomery multiplication**: Full implementation with software reduction
- **FJ64 hash table**: Uploaded to GPU memory (512KB)
- **Trial division**: 30 primes, same as CPU version

### Performance Notes
- GPU throughput: ~700K n/sec at n=10^12 (M4 Max)
- CPU throughput: ~3M n/sec single-threaded at n=10^12
- The software 128-bit arithmetic overhead (~5x) currently makes GPU slower than CPU
- GPU is still useful for verification or if CPU is busy with other tasks

## Analysis Tools

The `analysis/` directory contains profiling and tuning utilities:
- `benchmark_montgomery.c` - Compare Montgomery vs standard arithmetic
- `trial_div_tuning.c` - Optimize trial division prime count
- `profile_breakdown.c` - Time breakdown by function
