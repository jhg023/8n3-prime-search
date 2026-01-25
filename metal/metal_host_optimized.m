/*
 * Optimized Metal Host Implementation
 *
 * Key features:
 * - FJ64 hash table for deterministic 2-witness primality test
 * - 32-bit fast path for small primes (avoids 128-bit math)
 * - 100-prime trial division for better composite filtering
 */

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "metal_host.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mach/mach_time.h>

/* ========================================================================== */
/* Global State                                                                */
/* ========================================================================== */

static id<MTLDevice> g_device = nil;
static id<MTLCommandQueue> g_commandQueue = nil;
static id<MTLComputePipelineState> g_searchPipeline = nil;
static id<MTLComputePipelineState> g_searchPipelineOptimized = nil;
static id<MTLBuffer> g_fj64Buffer = nil;  /* FJ64 witness table */
static char g_deviceName[256] = {0};
static uint32_t g_maxBatchSize = 65536;
static bool g_useOptimized = true;  /* Use optimized kernel by default */

/* Statistics */
static GPUStats g_stats = {0};

/* Timing */
static mach_timebase_info_data_t g_timebaseInfo;

/* ========================================================================== */
/* Helper Functions                                                            */
/* ========================================================================== */

static double ticksToMilliseconds(uint64_t ticks) {
    if (g_timebaseInfo.denom == 0) {
        mach_timebase_info(&g_timebaseInfo);
    }
    return (double)ticks * g_timebaseInfo.numer / g_timebaseInfo.denom / 1000000.0;
}

/* ========================================================================== */
/* Initialization and Cleanup                                                  */
/* ========================================================================== */

bool metal_init(const uint16_t* fj64_bases) {
    @autoreleasepool {
        g_device = MTLCreateSystemDefaultDevice();
        if (!g_device) {
            fprintf(stderr, "Metal: Failed to get default device\n");
            return false;
        }

        strncpy(g_deviceName, [[g_device name] UTF8String], sizeof(g_deviceName) - 1);
        g_deviceName[sizeof(g_deviceName) - 1] = '\0';

        g_commandQueue = [g_device newCommandQueue];
        if (!g_commandQueue) {
            fprintf(stderr, "Metal: Failed to create command queue\n");
            return false;
        }

        /* Upload FJ64 witness table to GPU (512KB = 262144 * 2 bytes) */
        if (fj64_bases) {
            size_t fj64Size = 262144 * sizeof(uint16_t);
            g_fj64Buffer = [g_device newBufferWithBytes:fj64_bases
                                                 length:fj64Size
                                                options:MTLResourceStorageModeShared];
            if (!g_fj64Buffer) {
                fprintf(stderr, "Metal: Failed to upload FJ64 table\n");
                return false;
            }
        } else {
            fprintf(stderr, "Metal: Warning - no FJ64 table provided\n");
        }

        NSError* error = nil;
        id<MTLLibrary> library = nil;

        /* Try to load pre-compiled metallib first */
        NSString* libPath = @"metal/prime_search_optimized.metallib";
        if ([[NSFileManager defaultManager] fileExistsAtPath:libPath]) {
            NSURL* libURL = [NSURL fileURLWithPath:libPath];
            library = [g_device newLibraryWithURL:libURL error:&error];
        }

        if (!library) {
            /* Compile from source */
            NSString* sourcePath = @"metal/prime_search_optimized.metal";
            if ([[NSFileManager defaultManager] fileExistsAtPath:sourcePath]) {
                NSString* source = [NSString stringWithContentsOfFile:sourcePath
                                                             encoding:NSUTF8StringEncoding
                                                                error:&error];
                if (source) {
                    MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
                    if (@available(macOS 15.0, *)) {
                        options.mathMode = MTLMathModeFast;
                    } else {
                        #pragma clang diagnostic push
                        #pragma clang diagnostic ignored "-Wdeprecated-declarations"
                        options.fastMathEnabled = YES;
                        #pragma clang diagnostic pop
                    }
                    library = [g_device newLibraryWithSource:source options:options error:&error];
                }
            }
        }

        if (!library) {
            fprintf(stderr, "Metal: Failed to load shader library: %s\n",
                    error ? [[error localizedDescription] UTF8String] : "unknown error");
            return false;
        }

        /* Get optimized kernel */
        id<MTLFunction> optimizedFunc = [library newFunctionWithName:@"search_kernel_optimized"];
        if (optimizedFunc) {
            g_searchPipelineOptimized = [g_device newComputePipelineStateWithFunction:optimizedFunc error:&error];
            if (!g_searchPipelineOptimized) {
                fprintf(stderr, "Metal: Warning - optimized kernel failed: %s\n",
                        [[error localizedDescription] UTF8String]);
            }
        }

        /* Get standard kernel as fallback */
        id<MTLFunction> standardFunc = [library newFunctionWithName:@"search_kernel"];
        if (standardFunc) {
            g_searchPipeline = [g_device newComputePipelineStateWithFunction:standardFunc error:&error];
        }

        if (!g_searchPipelineOptimized && !g_searchPipeline) {
            fprintf(stderr, "Metal: Failed to create any pipeline state\n");
            return false;
        }

        /* Use optimized if available */
        g_useOptimized = (g_searchPipelineOptimized != nil);

        mach_timebase_info(&g_timebaseInfo);
        metal_reset_stats();

        return true;
    }
}

void metal_cleanup(void) {
    @autoreleasepool {
        g_fj64Buffer = nil;
        g_searchPipelineOptimized = nil;
        g_searchPipeline = nil;
        g_commandQueue = nil;
        g_device = nil;
        g_deviceName[0] = '\0';
    }
}

bool metal_is_available(void) {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return device != nil;
    }
}

const char* metal_get_device_name(void) {
    return g_deviceName;
}

/* ========================================================================== */
/* Search Operations                                                           */
/* ========================================================================== */

uint64_t metal_search_batch(
    const uint64_t* n_values,
    uint32_t count,
    GPUSearchResult* results,
    uint32_t threads_per_group
) {
    if (!g_device || count == 0) {
        return 0;
    }

    @autoreleasepool {
        /* Select pipeline */
        id<MTLComputePipelineState> pipeline = g_useOptimized ? g_searchPipelineOptimized : g_searchPipeline;
        if (!pipeline) {
            fprintf(stderr, "Metal: No pipeline available\n");
            return 0;
        }

        /* Create input buffer for n values */
        size_t nBufferSize = count * sizeof(uint64_t);
        id<MTLBuffer> nBuffer = [g_device newBufferWithBytes:n_values
                                                      length:nBufferSize
                                                     options:MTLResourceStorageModeShared];
        if (!nBuffer) {
            fprintf(stderr, "Metal: Failed to allocate n_values buffer\n");
            return 0;
        }

        /* Create output buffer for results */
        size_t resultBufferSize = count * sizeof(GPUSearchResult);
        id<MTLBuffer> resultBuffer = [g_device newBufferWithLength:resultBufferSize
                                                           options:MTLResourceStorageModeShared];
        if (!resultBuffer) {
            fprintf(stderr, "Metal: Failed to allocate results buffer\n");
            return 0;
        }

        /* Create batch size buffer */
        id<MTLBuffer> batchSizeBuffer = [g_device newBufferWithBytes:&count
                                                              length:sizeof(uint32_t)
                                                             options:MTLResourceStorageModeShared];

        /* Create command buffer */
        id<MTLCommandBuffer> commandBuffer = [g_commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:pipeline];
        [encoder setBuffer:nBuffer offset:0 atIndex:0];
        [encoder setBuffer:g_fj64Buffer offset:0 atIndex:1];  /* FJ64 witness table */
        [encoder setBuffer:resultBuffer offset:0 atIndex:2];
        [encoder setBuffer:batchSizeBuffer offset:0 atIndex:3];

        /* Configure thread groups */
        NSUInteger threadsPerThreadgroup = threads_per_group;
        NSUInteger maxThreadsPerGroup = [pipeline maxTotalThreadsPerThreadgroup];
        if (threadsPerThreadgroup > maxThreadsPerGroup) {
            threadsPerThreadgroup = maxThreadsPerGroup;
        }

        MTLSize threadgroupSize = MTLSizeMake(threadsPerThreadgroup, 1, 1);
        MTLSize gridSize = MTLSizeMake(threadsPerThreadgroup * count, 1, 1);

        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
        [encoder endEncoding];

        uint64_t startTime = mach_absolute_time();

        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        uint64_t endTime = mach_absolute_time();
        double elapsedMs = ticksToMilliseconds(endTime - startTime);

        if ([commandBuffer status] == MTLCommandBufferStatusError) {
            fprintf(stderr, "Metal: Command buffer error: %s\n",
                    [[[commandBuffer error] localizedDescription] UTF8String]);
            return 0;
        }

        memcpy(results, [resultBuffer contents], resultBufferSize);

        uint64_t counterexamples = 0;
        for (uint32_t i = 0; i < count; i++) {
            if (results[i].found == 0) {
                counterexamples++;
            }
        }

        g_stats.total_n_processed += count;
        g_stats.total_counterexamples += counterexamples;
        g_stats.total_gpu_time_ms += elapsedMs;
        g_stats.total_batches++;

        return counterexamples;
    }
}

/* ========================================================================== */
/* Statistics                                                                  */
/* ========================================================================== */

void metal_get_stats(GPUStats* stats) {
    if (stats) {
        *stats = g_stats;
    }
}

void metal_reset_stats(void) {
    memset(&g_stats, 0, sizeof(g_stats));
}

/* ========================================================================== */
/* Configuration                                                               */
/* ========================================================================== */

void metal_set_max_batch_size(uint32_t max_batch_size) {
    if (max_batch_size > 0) {
        g_maxBatchSize = max_batch_size;
    }
}

uint32_t metal_get_recommended_batch_size(void) {
    if (!g_device) {
        return 65536;
    }

    @autoreleasepool {
        NSUInteger recommendedMemory = [g_device recommendedMaxWorkingSetSize];
        uint32_t memoryBasedSize = (uint32_t)(recommendedMemory / 512);
        if (memoryBasedSize > 1000000) memoryBasedSize = 1000000;
        if (memoryBasedSize < 10000) memoryBasedSize = 10000;
        return memoryBasedSize;
    }
}

uint32_t metal_get_compute_units(void) {
    if (!g_device) {
        return 0;
    }

    @autoreleasepool {
        id<MTLComputePipelineState> pipeline = g_useOptimized ? g_searchPipelineOptimized : g_searchPipeline;
        if (!pipeline) return 0;

        NSUInteger maxThreadsPerGroup = [pipeline maxTotalThreadsPerThreadgroup];
        NSUInteger execWidth = [pipeline threadExecutionWidth];
        return (uint32_t)(maxThreadsPerGroup / execWidth);
    }
}
