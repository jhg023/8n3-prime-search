/*
 * Metal Host Implementation for GPU Prime Search
 *
 * Objective-C implementation of the Metal compute pipeline.
 * Handles device setup, buffer management, and kernel dispatch.
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
static id<MTLBuffer> g_fj64Buffer = nil;
static char g_deviceName[256] = {0};
static uint32_t g_maxBatchSize = 65536;  /* Default 64K */

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
        /* Get default Metal device */
        g_device = MTLCreateSystemDefaultDevice();
        if (!g_device) {
            fprintf(stderr, "Metal: Failed to get default device\n");
            return false;
        }

        /* Store device name */
        strncpy(g_deviceName, [[g_device name] UTF8String], sizeof(g_deviceName) - 1);
        g_deviceName[sizeof(g_deviceName) - 1] = '\0';

        /* Create command queue */
        g_commandQueue = [g_device newCommandQueue];
        if (!g_commandQueue) {
            fprintf(stderr, "Metal: Failed to create command queue\n");
            return false;
        }

        /* Load the Metal library */
        NSError* error = nil;
        id<MTLLibrary> library = nil;

        /* Try to load pre-compiled metallib first */
        NSString* libPath = [[NSBundle mainBundle] pathForResource:@"prime_search" ofType:@"metallib"];
        if (!libPath) {
            /* Try current directory */
            libPath = @"metal/prime_search.metallib";
        }

        if ([[NSFileManager defaultManager] fileExistsAtPath:libPath]) {
            NSURL* libURL = [NSURL fileURLWithPath:libPath];
            library = [g_device newLibraryWithURL:libURL error:&error];
        }

        if (!library) {
            /* Try loading from source */
            NSString* sourcePath = @"metal/prime_search.metal";
            if (![[NSFileManager defaultManager] fileExistsAtPath:sourcePath]) {
                /* Try bundle path */
                sourcePath = [[NSBundle mainBundle] pathForResource:@"prime_search" ofType:@"metal"];
            }

            if (sourcePath && [[NSFileManager defaultManager] fileExistsAtPath:sourcePath]) {
                NSString* source = [NSString stringWithContentsOfFile:sourcePath
                                                             encoding:NSUTF8StringEncoding
                                                                error:&error];
                if (source) {
                    MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
                    /* Use mathMode for macOS 15+ compatibility */
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

        /* Get the search kernel function */
        id<MTLFunction> searchFunction = [library newFunctionWithName:@"search_kernel"];
        if (!searchFunction) {
            fprintf(stderr, "Metal: Failed to find search_kernel function\n");
            return false;
        }

        /* Create compute pipeline state */
        g_searchPipeline = [g_device newComputePipelineStateWithFunction:searchFunction error:&error];
        if (!g_searchPipeline) {
            fprintf(stderr, "Metal: Failed to create pipeline state: %s\n",
                    [[error localizedDescription] UTF8String]);
            return false;
        }

        /* Upload fj64_bases table to GPU */
        size_t fj64Size = 262144 * sizeof(uint16_t);  /* 512KB */
        g_fj64Buffer = [g_device newBufferWithBytes:fj64_bases
                                             length:fj64Size
                                            options:MTLResourceStorageModeShared];
        if (!g_fj64Buffer) {
            fprintf(stderr, "Metal: Failed to allocate fj64_bases buffer\n");
            return false;
        }

        /* Initialize timing */
        mach_timebase_info(&g_timebaseInfo);

        /* Reset statistics */
        metal_reset_stats();

        return true;
    }
}

void metal_cleanup(void) {
    @autoreleasepool {
        g_fj64Buffer = nil;
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
    if (!g_device || !g_searchPipeline || !g_fj64Buffer || count == 0) {
        return 0;
    }

    @autoreleasepool {
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
        if (!batchSizeBuffer) {
            fprintf(stderr, "Metal: Failed to allocate batch_size buffer\n");
            return 0;
        }

        /* Create command buffer */
        id<MTLCommandBuffer> commandBuffer = [g_commandQueue commandBuffer];
        if (!commandBuffer) {
            fprintf(stderr, "Metal: Failed to create command buffer\n");
            return 0;
        }

        /* Create compute command encoder */
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        if (!encoder) {
            fprintf(stderr, "Metal: Failed to create compute encoder\n");
            return 0;
        }

        /* Set pipeline and buffers */
        [encoder setComputePipelineState:g_searchPipeline];
        [encoder setBuffer:nBuffer offset:0 atIndex:0];
        [encoder setBuffer:g_fj64Buffer offset:0 atIndex:1];
        [encoder setBuffer:resultBuffer offset:0 atIndex:2];
        [encoder setBuffer:batchSizeBuffer offset:0 atIndex:3];

        /* Configure thread groups */
        /* One threadgroup per n value, threads_per_group threads each */
        NSUInteger threadsPerThreadgroup = threads_per_group;
        NSUInteger maxThreadsPerGroup = [g_searchPipeline maxTotalThreadsPerThreadgroup];
        if (threadsPerThreadgroup > maxThreadsPerGroup) {
            threadsPerThreadgroup = maxThreadsPerGroup;
        }

        MTLSize threadgroupSize = MTLSizeMake(threadsPerThreadgroup, 1, 1);
        MTLSize gridSize = MTLSizeMake(threadsPerThreadgroup * count, 1, 1);

        /* Dispatch */
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
        [encoder endEncoding];

        /* Record start time */
        uint64_t startTime = mach_absolute_time();

        /* Submit and wait */
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        /* Record end time */
        uint64_t endTime = mach_absolute_time();
        double elapsedMs = ticksToMilliseconds(endTime - startTime);

        /* Check for errors */
        if ([commandBuffer status] == MTLCommandBufferStatusError) {
            fprintf(stderr, "Metal: Command buffer error: %s\n",
                    [[[commandBuffer error] localizedDescription] UTF8String]);
            return 0;
        }

        /* Copy results back */
        memcpy(results, [resultBuffer contents], resultBufferSize);

        /* Count counterexamples and update stats */
        uint64_t counterexamples = 0;
        for (uint32_t i = 0; i < count; i++) {
            if (results[i].found == 0) {
                counterexamples++;
            }
        }

        /* Update global statistics */
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
        return 65536;  /* Default */
    }

    @autoreleasepool {
        /* Get available memory (simplified heuristic) */
        /* On Apple Silicon, we can use more memory efficiently */
        NSUInteger recommendedMemory = [g_device recommendedMaxWorkingSetSize];

        /* Each n value uses:
         * - 8 bytes in n_values buffer
         * - 32 bytes in results buffer
         * - ~256 bytes threadgroup memory
         * Total: ~300 bytes per n value
         */
        uint32_t memoryBasedSize = (uint32_t)(recommendedMemory / 1024);  /* Leave headroom */
        if (memoryBasedSize > 1000000) memoryBasedSize = 1000000;  /* Cap at 1M */
        if (memoryBasedSize < 10000) memoryBasedSize = 10000;      /* Minimum 10K */

        return memoryBasedSize;
    }
}

uint32_t metal_get_compute_units(void) {
    if (!g_device) {
        return 0;
    }

    @autoreleasepool {
        /* Apple doesn't expose compute units directly, but we can estimate */
        /* Based on max threadgroup size and concurrent execution width */
        NSUInteger maxThreadsPerGroup = [g_searchPipeline maxTotalThreadsPerThreadgroup];
        NSUInteger execWidth = [g_searchPipeline threadExecutionWidth];

        /* Rough estimate: GPU cores ~ maxThreads / execWidth * some factor */
        /* This is imprecise but gives a useful relative measure */
        return (uint32_t)(maxThreadsPerGroup / execWidth);
    }
}
