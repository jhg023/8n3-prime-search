/*
 * Metal Host Interface for GPU Prime Search
 *
 * C-compatible interface for Metal GPU compute operations.
 * Provides initialization, batch search, and cleanup functions.
 */

#ifndef METAL_HOST_H
#define METAL_HOST_H

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ========================================================================== */
/* Types                                                                       */
/* ========================================================================== */

/* Result for a single n value */
typedef struct {
    uint64_t n;          /* The n value */
    uint64_t a;          /* Solution a (0 if counterexample) */
    uint64_t p;          /* Solution p (0 if counterexample) */
    uint32_t found;      /* 1 if solution found, 0 if counterexample */
    uint32_t _pad;       /* Padding for alignment */
} GPUSearchResult;

/* GPU statistics */
typedef struct {
    uint64_t total_n_processed;      /* Total n values processed */
    uint64_t total_counterexamples;  /* Number of potential counterexamples */
    double   total_gpu_time_ms;      /* Cumulative GPU execution time */
    uint64_t total_batches;          /* Number of batches executed */
} GPUStats;

/* ========================================================================== */
/* Initialization and Cleanup                                                  */
/* ========================================================================== */

/**
 * Initialize the Metal compute pipeline.
 *
 * This function:
 * - Creates the Metal device (GPU)
 * - Loads and compiles the compute shader
 * - Creates command queue and pipeline state
 * - Uploads the fj64_bases table to GPU memory
 *
 * @param fj64_bases  Pointer to the 512KB FJ64 witness table (262144 uint16_t values)
 * @return true on success, false on failure
 */
bool metal_init(const uint16_t* fj64_bases);

/**
 * Clean up Metal resources.
 *
 * Releases all GPU buffers, pipeline state, and device references.
 */
void metal_cleanup(void);

/**
 * Check if Metal GPU is available.
 *
 * @return true if a Metal-compatible GPU is available
 */
bool metal_is_available(void);

/**
 * Get GPU device name.
 *
 * @return Pointer to device name string (do not free)
 */
const char* metal_get_device_name(void);

/* ========================================================================== */
/* Search Operations                                                           */
/* ========================================================================== */

/**
 * Search a batch of n values on the GPU.
 *
 * For each n value, attempts to find (a, p) such that 8n + 3 = a^2 + 2p.
 * Uses the largest-a-first strategy (tests smallest primes first).
 *
 * @param n_values     Array of n values to search
 * @param count        Number of n values in the batch
 * @param results      Output array of GPUSearchResult (must be pre-allocated)
 * @param threads_per_group  Number of threads per threadgroup (typically 256)
 * @return Number of potential counterexamples found (results with found=0)
 */
uint64_t metal_search_batch(
    const uint64_t* n_values,
    uint32_t count,
    GPUSearchResult* results,
    uint32_t threads_per_group
);

/**
 * Get current GPU statistics.
 *
 * @param stats  Pointer to GPUStats structure to fill
 */
void metal_get_stats(GPUStats* stats);

/**
 * Reset GPU statistics counters.
 */
void metal_reset_stats(void);

/* ========================================================================== */
/* Configuration                                                               */
/* ========================================================================== */

/**
 * Set the maximum batch size for GPU operations.
 *
 * Larger batches have less overhead but use more GPU memory.
 * Default is 65536 (64K) n values per batch.
 *
 * @param max_batch_size  Maximum number of n values per GPU dispatch
 */
void metal_set_max_batch_size(uint32_t max_batch_size);


/**
 * Get the recommended batch size for the current GPU.
 *
 * @return Recommended batch size based on GPU memory and compute units
 */
uint32_t metal_get_recommended_batch_size(void);

/**
 * Get the number of GPU compute units.
 *
 * @return Number of compute units (execution units) on the GPU
 */
uint32_t metal_get_compute_units(void);

#ifdef __cplusplus
}
#endif

#endif /* METAL_HOST_H */
