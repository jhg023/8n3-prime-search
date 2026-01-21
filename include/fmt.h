/*
 * Formatting Utilities
 *
 * Number formatting functions for display output.
 */

#ifndef FMT_H
#define FMT_H

#include <stdint.h>

#define NUM_FMT_BUFS 8
#define NUM_FMT_SIZE 32

static char fmt_bufs[NUM_FMT_BUFS][NUM_FMT_SIZE];
static int fmt_buf_idx = 0;

/**
 * Format a number with comma separators (e.g., 1000000 -> "1,000,000")
 * Uses rotating static buffers - safe for up to NUM_FMT_BUFS concurrent uses
 */
static inline const char* fmt_num(uint64_t n) {
    char* buf = fmt_bufs[fmt_buf_idx];
    fmt_buf_idx = (fmt_buf_idx + 1) % NUM_FMT_BUFS;

    char temp[NUM_FMT_SIZE];
    int len = 0;
    if (n == 0) {
        temp[len++] = '0';
    } else {
        while (n > 0) {
            temp[len++] = '0' + (n % 10);
            n /= 10;
        }
    }

    int pos = 0;
    for (int i = len - 1; i >= 0; i--) {
        buf[pos++] = temp[i];
        if (i > 0 && i % 3 == 0) {
            buf[pos++] = ',';
        }
    }
    buf[pos] = '\0';

    return buf;
}

#endif /* FMT_H */
