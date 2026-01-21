# Makefile for Counterexample Search: 8n + 3 = a² + 2p

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -Iinclude
LDFLAGS = -lm

# Optimization flags
OPT_FLAGS = -O3 -march=native -mtune=native -flto -fomit-frame-pointer -funroll-loops -DNDEBUG

# Debug flags
DEBUG_FLAGS = -g -O0 -DDEBUG -fsanitize=address -fsanitize=undefined

# Source files
SRC_DIR = src
INCLUDE_DIR = include
BENCHMARK_DIR = benchmark

# Targets
TARGET = search
BENCHMARK_TARGET = benchmark

# Source files
SEARCH_SRC = $(SRC_DIR)/search.c
BENCHMARK_SRC = $(BENCHMARK_DIR)/benchmark.c

.PHONY: all release debug clean benchmark help

# Default: optimized build
all: release

# Optimized release build
release: CFLAGS += $(OPT_FLAGS)
release: $(TARGET)

# Debug build with sanitizers
debug: CFLAGS += $(DEBUG_FLAGS)
debug: LDFLAGS += -fsanitize=address -fsanitize=undefined
debug: $(TARGET)

# Main search program
$(TARGET): $(SEARCH_SRC) $(INCLUDE_DIR)/fj64_table.h
	$(CC) $(CFLAGS) -o $@ $(SEARCH_SRC) $(LDFLAGS)

# Benchmark program
benchmark: CFLAGS += $(OPT_FLAGS)
benchmark: $(BENCHMARK_SRC) $(INCLUDE_DIR)/fj64_table.h
	$(CC) $(CFLAGS) -o $(BENCHMARK_DIR)/$(BENCHMARK_TARGET) $(BENCHMARK_SRC) $(LDFLAGS)

# Clean build artifacts
clean:
	rm -f $(TARGET)
	rm -f $(BENCHMARK_DIR)/$(BENCHMARK_TARGET)
	rm -f *.o

# Run a quick test
test: release
	./$(TARGET) 1 10000

# Run benchmark
run-benchmark: benchmark
	./$(BENCHMARK_DIR)/$(BENCHMARK_TARGET)

# Help
help:
	@echo "Counterexample Search: 8n + 3 = a² + 2p"
	@echo ""
	@echo "Targets:"
	@echo "  all          Build optimized release (default)"
	@echo "  release      Build with full optimizations"
	@echo "  debug        Build with debug symbols and sanitizers"
	@echo "  benchmark    Build the benchmark program"
	@echo "  clean        Remove build artifacts"
	@echo "  test         Run a quick test (n = 1 to 10000)"
	@echo "  run-benchmark Run performance benchmarks"
	@echo "  help         Show this help message"
	@echo ""
	@echo "Usage:"
	@echo "  make              # Build optimized version"
	@echo "  ./search          # Run with default range [10^12, 10^12 + 10^6)"
	@echo "  ./search 1e9 2e9  # Run with custom range"
