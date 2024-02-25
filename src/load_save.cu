// #define _GNU_SOURCE
//
#include "wfc.cuh"
#include "bitfield.cuh"
#include "utils.cuh"

#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/// With a block side of 8, we have blocks of 8*8 := 64, which is the number of bits in an uint64_t.
static const uint8_t BLOCK_SIDE_U64 = 8;

static void
trim(char *restrict str)
{
    unsigned long start = 0, end = strlen(str) - 1;

    while (isspace(str[start])) {
        start++;
    }

    while (end > start && isspace(str[end])) {
        end--;
    }

    if (start > 0 || end < (strlen(str) - 1)) {
        memmove(str, str + start, end - start + 1);
        str[end - start + 1] = '\0';
    }
}

static char *
next(char *restrict str, char sep)
{
    char *ret = strchr(str, sep);
    if (NULL == ret) {
        fprintf(stderr, "failed to find character '%c'\n", sep);
        exit(EXIT_FAILURE);
    }
    ret[0] = '\0';
    ret += 1;
    return ret;
}

// static inline wfc_blocks *
// safe_malloc(uint64_t blkcnt)
// {
//     uint64_t size   = sizeof(wfc_blocks) + sizeof(uint64_t) * blkcnt;
//     wfc_blocks *ret = (wfc_blocks *)malloc(size);
//     if (ret != NULL) {
//         return ret;
//     } else {
//         fprintf(stderr, "failed to malloc %zu bytes\n", size);
//         exit(EXIT_FAILURE);
//     }
// }

static inline uint32_t
to_u32(const char *string)
{
    char *end          = NULL;
    const long integer = strtol(string, &end, 10);
    if (integer < 0) {
        fprintf(stderr, "expected positive integer, got %ld\n", integer);
        exit(EXIT_FAILURE);
    }
    return (uint32_t)integer;
}

static inline uint64_t
to_u64(const char *string)
{
    char *end               = NULL;
    const long long integer = strtoll(string, &end, 10);
    if (integer < 0) {
        fprintf(stderr, "expected positive integer, got %lld\n", integer);
        exit(EXIT_FAILURE);
    }
    return (uint64_t)integer;
}

wfc_blocks_ptr
wfc_load(uint64_t seed, const char *path)
{
    srandom((uint32_t)seed);

    ssize_t read    = -1;
    char *line      = NULL;
    size_t len      = 0;
    wfc_blocks *ret = NULL;
    uint64_t blkcnt = 0;

    FILE *restrict const f = fopen(path, "r");
    if (NULL == f) {
        fprintf(stderr, "failed to open `%s`: %s\n", path, strerror(errno));
        exit(EXIT_FAILURE);
    }

    if ((read = getline(&line, &len, f)) != -1) {
        const uint32_t block_side = to_u32(&line[1]);
        if (block_side > BLOCK_SIDE_U64) {
            fprintf(stderr, "invalid header of .dat file\n");
            exit(EXIT_FAILURE);
        }

        if (line[0] == 's') {
            blkcnt          = block_side * block_side;
            // ret             = safe_malloc(blkcnt + wfc_control_states_count(1, block_side));
            ret             = super_safe_malloc(1, block_side);
            ret->block_side = (uint8_t)block_side;
            ret->grid_side  = 1u;
            ret->seed = seed;
        } else if (line[0] == 'g') {
            blkcnt          = block_side * block_side;
            blkcnt          = blkcnt * blkcnt;
            // ret             = safe_malloc(blkcnt + wfc_control_states_count(block_side, block_side));
            ret             = super_safe_malloc(block_side, block_side);
            ret->block_side = (uint8_t)block_side;
            ret->grid_side  = (uint8_t)block_side;
            ret->seed = seed;
        } else {
            fprintf(stderr, "invalid header of .dat file\n");
            exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "invalid header of .dat file\n");
        exit(EXIT_FAILURE);
    }

    // TODO
    // {
    //     uint64_t mask       = 0;
    //     const uint64_t base = wfc_control_states_count(ret->grid_side, ret->block_side);
    //     for (uint8_t i = 0; i < ret->block_side * ret->block_side; i += 1) {
    //         mask = bitfield_set(mask, i);
    //     }
    //     for (uint64_t i = 0; i < blkcnt + base; i += 1) {
    //         ret->states[i] = mask;
    //     }
    // }

    while ((read = getline(&line, &len, f)) != -1) {
        trim(line);

        char *str_gy      = line;
        char *str_gx      = next(str_gy, ',');
        char *str_y       = next(str_gx, ',');
        char *str_x       = next(str_y, ',');
        char *str_state   = next(str_x, '=');
        const uint32_t gx = to_u32(str_gx), gy = to_u32(str_gy), x = to_u32(str_x),
                       y = to_u32(str_y);

        if (gx >= ret->grid_side || gy >= ret->grid_side) {
            fprintf(stderr, "invalid grid coordinates (%u, %u)\n", gx, gy);
            exit(EXIT_FAILURE);
        } else if (x >= ret->block_side || y >= ret->block_side) {
            fprintf(stderr, "invalid block coordinates (%u, %u) in grid (%u, %u)\n", x, y, gx, gy);
            exit(EXIT_FAILURE);
        }
      
        const uint64_t state   = to_u64(str_state);
        uint64_t collapsed = bitfield_set(0, state-1);
        // printf("collapsed (input): %lu (", state);
        // printBinary2(collapsed);
        // printf(") at : [%u, %u] [%u, %u]\n", gy, gx, y, x);
        host_grd_propagate_all(ret, gx, gy, x, y, collapsed);
        // grd_print(NULL, ret);
    }

    free(line);
    fclose(f);
    return ret;
}

void
wfc_save_into(const wfc_blocks_ptr blocks, const char data[], const char folder[])
{
    char destination[1024] = { 0 };
    const size_t data_len  = strlen(data);
    const char *file_name  = &data[data_len - 1];
    while (file_name != data && file_name[0] != '/') {
        file_name -= 1;
    }
    const char *file_end = strchr(file_name, '.');
    long length          = (file_end - file_name);
    if (length >= 1024) {
        length = 1023;
    } else if (length < 0) {
        length = 0;
    }

    const size_t folder_len = strlen(folder);
    if (folder[folder_len - 1] == '/' && file_name[0] == '/') {
        snprintf(destination, 1023, "%.*s%.*s.%lu.save", (int)(folder_len - 1), folder, (int)length,
                 file_name, blocks->seed);
    } else if ((folder[folder_len - 1] == '/' && file_name[0] != '/') ||
               (folder[folder_len - 1] != '/' && file_name[0] == '/')) {
        snprintf(destination, 1023, "%s%.*s.%lu.save", folder, (int)length, file_name,
                 blocks->seed);
    } else {
        snprintf(destination, 1023, "%s/%.*s.%lu.save", folder, (int)length, file_name,
                 blocks->seed);
    }
    fprintf(stdout, "save result to file: %s\n", destination);

    FILE *restrict f = fopen(destination, "w");
    if (NULL == f) {
        fprintf(stderr, "failed to open file: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    if (fprintf(f, "grid:  %hhu\n", blocks->grid_side) < 0) {
        fprintf(stderr, "failed to write: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }
    if (fprintf(f, "block: %hhu\n", blocks->block_side) < 0) {
        fprintf(stderr, "failed to write: %s\n", strerror(errno));
        exit(EXIT_FAILURE);
    }

    const uint64_t ends = blocks->grid_side * blocks->grid_side * blocks->block_side * blocks->block_side;

    for (uint64_t i = 0; i < ends; i += 1) {
        if (fprintf(f, "%lu ", (uint64_t)log2((double)blocks->states[i])+1) < 0) {
            fprintf(stderr, "failed to write: %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (i%(blocks->block_side * blocks->block_side) == (blocks->block_side * blocks->block_side) - 1) {
            if (fprintf(f, "\n") < 0) {
                fprintf(stderr, "failed to write: %s\n", strerror(errno));
                exit(EXIT_FAILURE);
            }
        }

    }

    fprintf(stdout, "saved successfully %lu states\n", ends);
    fclose(f);
}



void printBinary(uint64_t number) {
    // Determine the number of bits in u64_t
    int numBits = sizeof(uint64_t) * 8;

    // Loop through each bit in the number, starting from the most significant bit
    for (int i = numBits - 1; i >= 0; i--) {
        // Use a bitwise AND operation to check the value of the current bit
        if ((number & (1ULL << i)) != 0) {
            printf("1");
        } else {
            printf("0");
        }

        // Add space for better readability
        if (i % 8 == 0) {
            printf(" ");
        }
    }

    printf("\n");
}

