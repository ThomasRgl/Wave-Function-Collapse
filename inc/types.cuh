#pragma once

#include <inttypes.h>
#include <stdbool.h>

#define restrict  

/// Opaque type to store the seeds to try for the solving process. You may push to it and pop from
/// it. You may not try to index it manually or free this structure, it will be automatically freed
/// when no more items are present inside it.
typedef struct seeds_list seeds_list;




typedef struct {
    uint32_t x, y;
} vec2;

typedef struct {
    vec2 location;
    uint8_t entropy;

    uint8_t _1;
    uint16_t _2;
} entropy_location;

// typedef struct {
//     uint64_t seed;
//     uint8_t block_side;
//     uint8_t grid_side;
//
//     uint8_t _1;
//     uint8_t _2;
//     uint32_t _3;
//
//     // uint64_t states[];
//     uint64_t *states;
// } wfc_blocks;

typedef struct {
    uint32_t gx;
    uint32_t gy;
    uint32_t x;
    uint32_t y;
} vec4;

typedef struct {
    uint64_t seed;
    uint8_t block_side;
    uint8_t grid_side;
    bool solved;

    uint64_t * row_masks;
    uint64_t * col_masks;
    uint64_t * blk_masks;

    vec4 * stack_cells;
    uint32_t stack_size;

    uint64_t *states;
} wfc_blocks;

typedef wfc_blocks *wfc_blocks_ptr;

typedef struct {
    const char *const data_file;
    const char *const output_folder;
    seeds_list *restrict seeds;
    const uint64_t parallel;
    bool (*const solver)(wfc_blocks_ptr);
} wfc_args;

typedef struct {
    const char *const name;
    bool (*function)(wfc_blocks_ptr);
} wfc_solver;
