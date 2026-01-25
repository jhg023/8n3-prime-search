# Makefile.metal - Metal GPU build rules for 8n+3 Prime Search
#
# This Makefile handles compilation of Metal shaders and linking
# with Objective-C host code.
#
# Include from main Makefile: include metal/Makefile.metal

# Metal compiler (xcrun finds the correct SDK path)
METAL = xcrun -sdk macosx metal
METALLIB = xcrun -sdk macosx metallib

# Metal compilation flags
METAL_FLAGS = -std=metal3.0 -O3 -ffast-math

# Directories
METAL_DIR = metal
METAL_SRC = $(METAL_DIR)/prime_search.metal
METAL_AIR = $(METAL_DIR)/prime_search.air
METAL_LIB = $(METAL_DIR)/prime_search.metallib
METAL_HOST = $(METAL_DIR)/metal_host.m
METAL_HEADER = $(METAL_DIR)/metal_host.h

# GPU search program
GPU_TARGET = search_gpu
GPU_SRC = src/search_gpu.c

# Frameworks for Metal/Objective-C
METAL_FRAMEWORKS = -framework Metal -framework Foundation
METAL_OBJ_FLAGS = -fobjc-arc

# Headers
HEADERS = include/fmt.h include/arith.h include/prime.h \
          include/solve.h include/fj64_table.h include/arith_montgomery.h \
          $(METAL_HEADER)

.PHONY: metal metal-shader clean-metal

# Build GPU search program
metal: CFLAGS += -O3 -march=native -DNDEBUG
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

# Clean Metal build artifacts
clean-metal:
	rm -f $(METAL_AIR) $(METAL_LIB) $(GPU_TARGET)
