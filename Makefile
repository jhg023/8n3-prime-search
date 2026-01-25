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
    # Metal GPU support is available on macOS
    METAL_AVAILABLE = 1
else
    # Linux: Standard OpenMP
    OPENMP_CFLAGS = -fopenmp
    OPENMP_LDFLAGS = -fopenmp
    # Metal not available on Linux
    METAL_AVAILABLE = 0
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
          $(INCLUDE_DIR)/solve.h $(INCLUDE_DIR)/fj64_table.h $(INCLUDE_DIR)/arith_montgomery.h \
          $(INCLUDE_DIR)/prime_sieve.h $(INCLUDE_DIR)/batch_sieve.h $(INCLUDE_DIR)/residue_analysis.h

# Additional source files
SEARCH_BATCHED_SRC = $(SRC_DIR)/search_batched.c
BENCHMARK_APPROACHES_SRC = $(BENCHMARK_DIR)/benchmark_approaches.c

.PHONY: all release debug clean benchmark test run-benchmark help single-threaded metal metal-shader clean-metal test-gpu search_batched benchmark-approaches

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

# Batched search program (single-threaded)
search_batched: CFLAGS += $(OPT_FLAGS)
search_batched: $(SEARCH_BATCHED_SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o search_batched $(SEARCH_BATCHED_SRC) $(LDFLAGS)

# Benchmark approaches comparison tool
benchmark-approaches: CFLAGS += $(OPT_FLAGS)
benchmark-approaches: $(BENCHMARK_APPROACHES_SRC) $(HEADERS)
	$(CC) $(CFLAGS) -o $(BENCHMARK_DIR)/benchmark_approaches $(BENCHMARK_APPROACHES_SRC) $(LDFLAGS)

# Run approaches benchmark
run-benchmark-approaches: benchmark-approaches
	./$(BENCHMARK_DIR)/benchmark_approaches

# Run quick approaches benchmark
run-benchmark-approaches-quick: benchmark-approaches
	./$(BENCHMARK_DIR)/benchmark_approaches --quick

# Clean build artifacts
clean: clean-metal
	rm -f $(TARGET)
	rm -f search_batched
	rm -f $(BENCHMARK_DIR)/$(BENCHMARK_TARGET)
	rm -f $(BENCHMARK_DIR)/benchmark_approaches
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
	@echo "  search_batched    Build batched search (segmented sieve)"
	@echo "  benchmark-approaches  Build optimization comparison benchmark"
	@echo "  metal             Build GPU-accelerated version (macOS only)"
	@echo "  clean             Remove build artifacts"
	@echo "  test              Run a quick test (n = 1 to 10000)"
	@echo "  test-gpu          Test GPU against CPU (macOS only)"
	@echo "  run-benchmark     Run benchmark (10M iterations/scale)"
	@echo "  run-benchmark-quick  Run quick benchmark (1M iterations/scale)"
	@echo "  run-benchmark-approaches  Run optimization comparison benchmark"
	@echo "  help              Show this help message"
	@echo ""
	@echo "CPU Usage:"
	@echo "  make                      # Build optimized parallel version"
	@echo "  ./search                  # Run with default range, all cores"
	@echo "  ./search 1e9 2e9          # Run with custom range"
	@echo "  ./search 1e9 2e9 --threads 4  # Run with 4 threads"
	@echo "  ./search 1e9 2e9 --sieve-threshold 1e8  # Use prime sieve (12MB)"
	@echo ""
	@echo "Batched Search:"
	@echo "  make search_batched       # Build batched search"
	@echo "  ./search_batched          # Run with default range"
	@echo "  ./search_batched 1e9 2e9 --batch-size 131072"
	@echo ""
	@echo "GPU Usage (macOS with Metal):"
	@echo "  make metal                # Build GPU-accelerated version"
	@echo "  ./search_gpu              # Run with default range on GPU"
	@echo "  ./search_gpu 1e9 2e9      # Run custom range on GPU"
	@echo "  ./search_gpu 1e9 2e9 --batch-size 100000"
	@echo ""
	@echo "Benchmark Options:"
	@echo "  ./benchmark/benchmark_suite             # Default (10M iterations)"
	@echo "  ./benchmark/benchmark_suite --quick     # Quick (1M iterations)"
	@echo "  ./benchmark/benchmark_suite --count N   # Custom iteration count"
	@echo "  ./benchmark/benchmark_approaches        # Compare optimization approaches"
	@echo ""
	@echo "Note: On macOS, install libomp via: brew install libomp"

# ============================================================================
# Metal GPU Support (macOS only)
# ============================================================================

ifeq ($(METAL_AVAILABLE),1)

# Metal compiler (xcrun finds the correct SDK path)
METAL = xcrun -sdk macosx metal
METALLIB = xcrun -sdk macosx metallib

# Metal compilation flags
METAL_FLAGS = -std=metal3.0 -O3 -ffast-math

# Metal directories and files
METAL_DIR = metal
METAL_SRC = $(METAL_DIR)/prime_search_optimized.metal
METAL_AIR = $(METAL_DIR)/prime_search_optimized.air
METAL_LIB = $(METAL_DIR)/prime_search_optimized.metallib
METAL_HOST = $(METAL_DIR)/metal_host_optimized.m
METAL_HEADER = $(METAL_DIR)/metal_host.h

# GPU search program
GPU_TARGET = search_gpu
GPU_SRC = src/search_gpu.c

# Frameworks for Metal/Objective-C
METAL_FRAMEWORKS = -framework Metal -framework Foundation
METAL_OBJ_FLAGS = -fobjc-arc

# Build GPU search program with Metal support
metal: CFLAGS += $(OPT_FLAGS)
metal: $(METAL_LIB) $(GPU_TARGET)

# Compile Metal shader: .metal -> .air
$(METAL_AIR): $(METAL_SRC)
	@echo "Compiling Metal shader..."
	$(METAL) $(METAL_FLAGS) -c $< -o $@

# Link Metal library: .air -> .metallib
$(METAL_LIB): $(METAL_AIR)
	@echo "Creating Metal library..."
	$(METALLIB) -o $@ $<

# Build GPU search program
$(GPU_TARGET): $(GPU_SRC) $(METAL_HOST) $(METAL_LIB) $(HEADERS)
	@echo "Building GPU search program..."
	$(CC) $(CFLAGS) -Iinclude -I$(METAL_DIR) \
		$(METAL_OBJ_FLAGS) \
		-o $@ $(GPU_SRC) $(METAL_HOST) \
		$(METAL_FRAMEWORKS) $(LDFLAGS)

# Compile Metal shader only (for testing)
metal-shader: $(METAL_LIB)

# Test GPU against CPU
test-gpu: metal
	./$(GPU_TARGET) --verify-only

# Clean Metal build artifacts
clean-metal:
	rm -f $(METAL_AIR) $(METAL_LIB) $(GPU_TARGET)

else
# Metal not available (non-macOS)
metal:
	@echo "Error: Metal GPU support is only available on macOS"
	@exit 1

metal-shader:
	@echo "Error: Metal GPU support is only available on macOS"
	@exit 1

test-gpu:
	@echo "Error: Metal GPU support is only available on macOS"
	@exit 1

clean-metal:
	@echo "Metal not available - nothing to clean"

endif
