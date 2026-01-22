# Makefile for Counterexample Search: 8n + 3 = a^2 + 2p

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -Iinclude
LDFLAGS = -lm

# OpenMP flags (auto-detect platform)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    # macOS: Use Homebrew's libomp if available
    BREW_PREFIX := $(shell brew --prefix 2>/dev/null || echo /usr/local)
    OPENMP_CFLAGS = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/opt/libomp/include
    OPENMP_LDFLAGS = -L$(BREW_PREFIX)/opt/libomp/lib -lomp
else
    # Linux: Standard OpenMP
    OPENMP_CFLAGS = -fopenmp
    OPENMP_LDFLAGS = -fopenmp
endif

# Optimization flags
OPT_FLAGS = -O3 -march=native -mtune=native -flto -fomit-frame-pointer -funroll-loops -DNDEBUG

# Debug flags
DEBUG_FLAGS = -g -O0 -DDEBUG -fsanitize=address -fsanitize=undefined

# Directories
SRC_DIR = src
INCLUDE_DIR = include
BENCHMARK_DIR = benchmark

# Targets
TARGET = search
BENCHMARK_TARGET = benchmark_suite

# Source files
SEARCH_SRC = $(SRC_DIR)/search.c
BENCHMARK_SRC = $(BENCHMARK_DIR)/benchmark_suite.c

# Header dependencies
HEADERS = $(INCLUDE_DIR)/fmt.h $(INCLUDE_DIR)/arith.h $(INCLUDE_DIR)/prime.h \
          $(INCLUDE_DIR)/solve.h $(INCLUDE_DIR)/fj64_table.h $(INCLUDE_DIR)/arith_montgomery.h

.PHONY: all release debug clean benchmark test run-benchmark help single-threaded

# Default: optimized parallel build
all: release

# Optimized release build (with OpenMP)
release: CFLAGS += $(OPT_FLAGS) $(OPENMP_CFLAGS)
release: LDFLAGS += $(OPENMP_LDFLAGS)
release: $(TARGET)

# Single-threaded build (no OpenMP)
single-threaded: CFLAGS += $(OPT_FLAGS)
single-threaded: $(TARGET)

# Debug build with sanitizers (no OpenMP - conflicts with sanitizers)
debug: CFLAGS += $(DEBUG_FLAGS)
debug: LDFLAGS += -fsanitize=address -fsanitize=undefined
debug: $(TARGET)

# Main search program
$(TARGET): $(SEARCH_SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $(SEARCH_SRC) $(LDFLAGS)

# Benchmark suite (single-threaded for consistent comparisons)
benchmark: CFLAGS += $(OPT_FLAGS)
benchmark: $(BENCHMARK_SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o $(BENCHMARK_DIR)/$(BENCHMARK_TARGET) $(BENCHMARK_SRC) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(TARGET)
	rm -f $(BENCHMARK_DIR)/$(BENCHMARK_TARGET)
	rm -f *.o

# Run a quick test
test: release
	./$(TARGET) 1 10000

# Run benchmark suite
run-benchmark: benchmark
	./$(BENCHMARK_DIR)/$(BENCHMARK_TARGET)

# Run quick benchmark
run-benchmark-quick: benchmark
	./$(BENCHMARK_DIR)/$(BENCHMARK_TARGET) --quick

# Help
help:
	@echo "Counterexample Search: 8n + 3 = a^2 + 2p"
	@echo ""
	@echo "Targets:"
	@echo "  all               Build optimized parallel release (default)"
	@echo "  release           Build with full optimizations + OpenMP"
	@echo "  single-threaded   Build optimized without OpenMP"
	@echo "  debug             Build with debug symbols and sanitizers"
	@echo "  benchmark         Build the benchmark suite"
	@echo "  clean             Remove build artifacts"
	@echo "  test              Run a quick test (n = 1 to 10000)"
	@echo "  run-benchmark     Run benchmark (10M iterations/scale)"
	@echo "  run-benchmark-quick  Run quick benchmark (1M iterations/scale)"
	@echo "  help              Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make                      # Build optimized parallel version"
	@echo "  ./search                  # Run with default range, all cores"
	@echo "  ./search 1e9 2e9          # Run with custom range"
	@echo "  ./search 1e9 2e9 --threads 4  # Run with 4 threads"
	@echo ""
	@echo "Benchmark Options:"
	@echo "  ./benchmark/benchmark_suite             # Default (10M iterations)"
	@echo "  ./benchmark/benchmark_suite --quick     # Quick (1M iterations)"
	@echo "  ./benchmark/benchmark_suite --count N   # Custom iteration count"
	@echo ""
	@echo "Note: On macOS, install libomp via: brew install libomp"
