# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
