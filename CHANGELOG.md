# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.3.2] - 2026-01-22

### Changed
- **Branchless exponentiation** in Miller-Rabin witness test (`arith_montgomery.h`)
  - Use conditional select instead of conditional branch in exponentiation loop
  - Avoids branch misprediction when exponent bits are unpredictable
  - ~3% speedup for Miller-Rabin tests

### Performance
- ~2,796,000 n/sec at n = 10^12 (+2.6% from v1.3.1)
- ~2,126,000 n/sec at n = 10^15 (+2.6% from v1.3.1)
- ~1,598,000 n/sec at n = 2×10^18 (+4.2% from v1.3.1)

### Not Implemented (analyzed but no benefit)
- **Branchless Montgomery reduction**: Branch version is 7% faster (branch is predictable)
- **FJ64 table prefetch**: Table fits in cache; prefetch adds overhead
- **Alternative compiler flags**: `-Ofast`, `-ffast-math`, etc. provide no improvement
- **Loop alignment attributes**: No measurable benefit

## [1.3.1] - 2026-01-22

### Changed
- **Montgomery constant caching** in `is_prime_fj64_fast()` (`prime.h`)
  - Pre-compute `n_inv` and `r_sq` once and reuse for both Miller-Rabin witnesses
  - Saves ~8-12 ns per candidate by avoiding redundant computation
  - Added `mr_witness_montgomery_cached()` function in `arith_montgomery.h`
- **Incremental candidate tracking with delta** in `find_solution()` (`solve.h`)
  - Track delta separately instead of computing `2*(a-1)` each iteration
  - Update loop: `candidate += delta; delta -= 4; a -= 2`
  - Reduces iteration overhead from ~1 ns/iter to ~0.3 ns/iter

### Performance
- Modest improvements (0-5%) depending on scale
- ~2,725,000 n/sec at n = 10^12
- ~1,535,000 n/sec at n = 2×10^18

### Not Implemented (analyzed but no benefit)
- **Wheel factorization**: Lookup overhead (~0.6 ns) negates benefit of skipping 53% of candidates
- **Fewer trial primes**: 30 remains optimal; fewer primes = more expensive MR tests
- **Loop unrolling**: Compiler already does an excellent job
- **Branch hints**: `__builtin_expect` provided no measurable improvement

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
  - `trial_div_tuning.c` - Trial division count optimization
  - `trial_div_final.c` - Final 30 vs 120 prime comparison

### Changed
- **prime.h** now uses Montgomery-accelerated Miller-Rabin by default
- **Reduced trial division from 120 to 30 primes** (3-127)
  - With Montgomery-accelerated MR, fewer trial primes is optimal
  - 30 primes filter ~80% of composites with minimal overhead
  - 120 primes filter ~85% but 90 extra modulo ops cost more than occasional MR
  - Benchmarked: 30 primes is 7-11% faster than 120 at large n
- Header dependency: `arith.h` → `arith_montgomery.h` → `prime.h` → `solve.h`

### Performance
- **2.56x total speedup** at 64-bit scale (n ~ 2×10^18)
- ~5,200,000 n/sec at n = 10^6 (was ~3,000,000)
- ~2,600,000 n/sec at n = 10^12 (was ~1,400,000)
- ~2,000,000 n/sec at n = 10^15 (was ~1,000,000)
- ~1,510,000 n/sec at n = 2×10^18 (was ~570,000)

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
- **Extended trial division in prime.h**: Updated from 30 to 120 primes (later reverted in v1.3.0)

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
