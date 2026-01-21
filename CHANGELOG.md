# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [1.1.0] - 2026-01-21

### Changed
- **Reverse iteration order**: Now tests largest `a` values first (smallest prime candidates)
  - Smaller primes are faster to verify and more likely to be prime
  - Improves solution-finding speed
- **Expanded trial division**: Increased from 30 primes (up to 127) to 120 primes (up to 661)
  - Benchmarked across different counts to find optimal configuration
  - +12.4% throughput improvement at large n (~2×10¹⁸)
- **Unbuffered stdout**: Added `setbuf(stdout, NULL)` for real-time progress output

### Added
- **32-bit candidate tracking**: Progress reports now show percentage of candidates fitting in 32-bit
- **Candidate statistics**: Final results include total candidates tested and 32-bit breakdown

### Fixed
- **Loop underflow bug**: Fixed unsigned integer underflow when iterating `a` values in reverse
  - Previous code could cause infinite loops or false counterexamples

### Performance
- ~1,360,000 n/sec throughput at n = 10¹² (100% 32-bit candidates)
- ~970,000 n/sec throughput at n = 10¹⁵ (76% 32-bit candidates)
- ~560,000 n/sec throughput at n = 2×10¹⁸ (6% 32-bit candidates)

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
- ~460,000 n/sec throughput at n = 10¹²
- ~380,000 n/sec throughput at n = 10¹⁴
- ~310,000 n/sec throughput at n = 10¹⁶

### Technical Details
- FJ64_262K uses 512KB hash table for witness selection
- Only 2 Miller-Rabin tests per candidate (down from 7)
- Trial division filters ~82% of composites before primality testing
