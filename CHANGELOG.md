# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.3.0] - 2026-01-21

### Added
- **Montgomery multiplication** for Miller-Rabin tests (`arith_montgomery.h`)
  - 3x faster modular arithmetic for moduli < 2^63
  - Hybrid implementation falls back to `__uint128_t` for n ≥ 2^63
  - Avoids expensive division instructions, uses multiplication + shifts
- **Analysis tools** for optimization research:
  - `profile_breakdown.c` - Time breakdown by component
  - `trial_div_profile.c` - Trial division hit rate analysis
  - `wheel_analysis.c` - Wheel factorization potential (54.5% skippable)
  - `benchmark_montgomery.c` - Montgomery vs standard comparison

### Changed
- **prime.h** now uses Montgomery-accelerated Miller-Rabin by default
- Header dependency: `arith.h` → `arith_montgomery.h` → `prime.h` → `solve.h`

### Performance
- **2.3x speedup** at 64-bit scale (n ~ 2×10^18)
- ~4,500,000 n/sec at n = 10^6 (was ~3,000,000)
- ~2,300,000 n/sec at n = 10^12 (was ~1,400,000)
- ~1,800,000 n/sec at n = 10^15 (was ~1,000,000)
- ~1,360,000 n/sec at n = 2×10^18 (was ~570,000)

## [1.2.0] - 2026-01-21

### Changed
- **Refactored codebase**: Deduplicated code by moving shared functionality to header files
  - `search.c` now uses shared headers instead of inline implementations
  - Reduced `search.c` from ~520 lines to ~290 lines
- **Merged solve headers**: Combined `solve.h` and `solve_opt.h` into single `solve.h`
  - All iteration strategies now in one file
  - Default `find_solution()` uses large-to-small (fastest strategy)
- **Simplified benchmark suite**: Replaced complex multi-mode benchmark with focused throughput test
  - Tests 7 scales from 10^6 to 2*10^18 (~64-bit limit)
  - Reports rate, average checks per n, and time
  - 10M iterations per scale for stable measurements
- **Updated default search range**: Changed from 1M to 10M iterations (10^12 to 10^12 + 10^7)
- **Extended trial division in prime.h**: Updated from 30 to 120 primes to match production code

### Added
- **Optional stats tracking in solve.h**: Define `SOLVE_TRACK_STATS` to enable candidate statistics
- **Helper functions in solve.h**: `trial_division_check()` and `is_candidate_prime()`

### Removed
- `include/solve_opt.h` (merged into `solve.h`)
- `benchmark/benchmark.c` (replaced by `benchmark_suite.c`)
- Duplicate implementations in `search.c`

### Fixed
- Header dependency chain now properly organized: `arith.h` -> `prime.h` -> `solve.h`

## [1.1.0] - 2026-01-21

### Changed
- **Reverse iteration order**: Now tests largest `a` values first (smallest prime candidates)
  - Smaller primes are faster to verify and more likely to be prime
  - Improves solution-finding speed
- **Expanded trial division**: Increased from 30 primes (up to 127) to 120 primes (up to 661)
  - Benchmarked across different counts to find optimal configuration
  - +12.4% throughput improvement at large n (~2*10^18)
- **Unbuffered stdout**: Added `setbuf(stdout, NULL)` for real-time progress output

### Added
- **32-bit candidate tracking**: Progress reports now show percentage of candidates fitting in 32-bit
- **Candidate statistics**: Final results include total candidates tested and 32-bit breakdown

### Fixed
- **Loop underflow bug**: Fixed unsigned integer underflow when iterating `a` values in reverse
  - Previous code could cause infinite loops or false counterexamples

### Performance
- ~1,360,000 n/sec throughput at n = 10^12 (100% 32-bit candidates)
- ~970,000 n/sec throughput at n = 10^15 (76% 32-bit candidates)
- ~560,000 n/sec throughput at n = 2*10^18 (6% 32-bit candidates)

## [1.0.0] - 2025-01-21

### Added
- Initial release
- FJ64_262K primality test integration (1.9x faster than 7-witness MR)
- Trial division pre-filtering by 30 small primes
- Support for scientific notation in command-line arguments
- Progress reporting with throughput statistics
- Verification of known solutions on startup
- Comprehensive benchmark suite comparing FJ64_262K, BPSW, and MR-7
- Documentation of algorithm and performance characteristics

### Performance
- ~460,000 n/sec throughput at n = 10^12
- ~380,000 n/sec throughput at n = 10^14
- ~310,000 n/sec throughput at n = 10^16

### Technical Details
- FJ64_262K uses 512KB hash table for witness selection
- Only 2 Miller-Rabin tests per candidate (down from 7)
- Trial division filters ~82% of composites before primality testing
