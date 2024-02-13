#define _GNU_SOURCE

#include "wfc.h"
#include "bitfield.h"
// #include "utils.h" 
#include "md5.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strings.h>

uint64_t
entropy_collapse_state(uint64_t state,
                       uint32_t gx, uint32_t gy, uint32_t x, uint32_t y,
                       uint64_t seed,
                       uint64_t iteration)
{
    uint8_t digest[16]     = { 0 };
    uint64_t random_number = 0;
    struct {
        uint32_t gx, gy, x, y;
        uint64_t seed, iteration;
    } random_state = {
        .gx        = gx,
        .gy        = gy,
        .x         = x,
        .y         = y,
        .seed      = seed,
        .iteration = iteration,
    };

    md5((uint8_t *)&random_state, sizeof(random_state), digest);

    return 0;
}

uint8_t
entropy_compute(uint64_t state)
{
    return 0;
}

void
wfc_clone_into(wfc_blocks_ptr *const restrict ret_ptr, uint64_t seed, const wfc_blocks_ptr blocks)
{
    const uint64_t grid_size  = blocks->grid_side;
    const uint64_t block_size = blocks->block_side;
    wfc_blocks_ptr ret        = *ret_ptr;

    const uint64_t size = (wfc_control_states_count(grid_size, block_size) * sizeof(uint64_t)) +
                          (grid_size * grid_size * block_size * block_size * sizeof(uint64_t)) +
                          sizeof(wfc_blocks);

    if (NULL == ret) {
        if (NULL == (ret = malloc(size))) {
            fprintf(stderr, "failed to clone blocks structure\n");
            exit(EXIT_FAILURE);
        }
    } else if (grid_size != ret->grid_side || block_size != ret->block_side) {
        fprintf(stderr, "size mismatch!\n");
        exit(EXIT_FAILURE);
    }

    memcpy(ret, blocks, size);
    ret->seed = seed;
    *ret_ptr       = ret;
}

entropy_location
blk_min_entropy(const wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy)
{
    vec2 the_location   = { 0 };
    uint8_t min_entropy = UINT8_MAX;

    entropy_location n;
    return n;
}

static inline uint64_t
blk_filter_mask_for_row(wfc_blocks_ptr blocks,
                           uint32_t gx, uint32_t gy, uint32_t y)
{
    uint64_t mask = 0;
    uint64_t state = 0;
    for (uint32_t i = 0; i < blocks->block_side; i++) {
        state = *blk_at(blocks, gx, gy, i, y);
        if(bitfield_count(state) == 1){
            mask |= state;
        }
    }

    return mask;
}

static inline uint64_t
blk_filter_mask_for_column(wfc_blocks_ptr blocks,
                            uint32_t gx, uint32_t gy, uint32_t x)
{
    uint64_t mask = 0;
    uint64_t state = 0;
    for (uint32_t i = 0; i < blocks->block_side; i++) {
        state = *blk_at(blocks, gx, gy, x, i);
        if(bitfield_count(state) == 1){
            mask |= state;
        }
    }

    return mask;
}

static inline uint64_t
grd_filter_mask_for_row(wfc_blocks_ptr blocks,
                           uint32_t gy, uint32_t y)
{
    uint64_t mask = 0;
    for (uint32_t i = 0; i < blocks->grid_side; i++) {
        mask |= blk_filter_mask_for_row(blocks, i, gy, y);
    }
    
    return mask;
}

static inline uint64_t
grd_filter_mask_for_column(wfc_blocks_ptr blocks,
                        uint32_t gx, uint32_t x)
{
    uint64_t mask = 0;
    for (uint32_t i = 0; i < blocks->grid_side; i++) {
        mask |= blk_filter_mask_for_column(blocks, gx, i, x);
    }

    return mask;
}

static inline uint64_t
blk_filter_mask_for_block(wfc_blocks_ptr blocks,
                          uint32_t gy, uint32_t gx)
{
    uint64_t mask = 0;
    uint64_t state = 0;
    for (uint32_t i = 0; i < blocks->block_side; i++) {
        for (uint32_t j = 0; j < blocks->block_side; j++) {
            state = *blk_at(blocks, gx, gy, i, j);
            if(bitfield_count(state) == 1){
                mask |= state;
            }
        }
    }

    return mask;
}

bool
grd_check_error_in_column(wfc_blocks_ptr blocks, uint32_t gx, u_int32_t gy, uint32_t x, u_int32_t y)
{
    uint64_t mask = grd_filter_mask_for_column(blocks, gx, x);
    uint64_t state = 0;
    for (uint32_t i = 0; i < blocks->grid_side; i++) {
        for (uint32_t j = 0; j < blocks->block_side; j++) {
            state = *blk_at(blocks, gx, i, x, j);
            if (bitfield_count(state) == 1) {
                if (bitfield_count(state|mask) != bitfield_count(mask)+1) {
                    // printf("Col : i: %u, j: %u, x: %u, y: %u, gx: %u, gy: %u\n", i, j, x, y, gx, gy);
                    if(i != gy && j != y) {
                        return true;
                    }
                }   
            }
        }
    }
  
    return false;
}

bool
grd_check_error_in_row(wfc_blocks_ptr blocks, uint32_t gx, u_int32_t gy, uint32_t x, u_int32_t y)
{
    uint64_t mask = grd_filter_mask_for_row(blocks, gy, y);
    uint64_t state = 0;
    for (uint32_t i = 0; i < blocks->grid_side; i++) {
        for (uint32_t j = 0; j < blocks->block_side; j++) {
            state = *blk_at(blocks, i, gy, j, y);
            if (bitfield_count(state) == 1) {
                if (bitfield_count(state|mask) != bitfield_count(mask)+1) {
                    // printf("Row : i: %u, j: %u, x: %u, y: %u, gx: %u, gy: %u\n", i, j, x, y, gx, gy);
                    if(i != gx && j != x) {
                        return true;
                    }
                }   
            }
        }
    }

    return false;
}

bool
grd_check_error_in_block(wfc_blocks_ptr blocks, uint32_t gx, u_int32_t gy, uint32_t x, u_int32_t y)
{
    uint64_t mask = blk_filter_mask_for_block(blocks, gx, gy);
    uint64_t state = 0;
    for (uint32_t i = 0; i < blocks->grid_side; i++) {
        for (uint32_t j = 0; j < blocks->grid_side; j++) {
            for (uint32_t bi = 0; bi < blocks->block_side; bi++) {
                for (uint32_t bj = 0; bj < blocks->block_side; bj++) {
                    state = *blk_at(blocks, i, j, bi, bj);
                    if (bitfield_count(state) == 1) {
                        if (bitfield_count(state|mask) != bitfield_count(mask)+1) {
                            // printf("Block : i: %u, j: %u, bi: %u, bj: %u, x: %u, y: %u, gx: %u, gy: %u\n", i, j, bi, bj, x, y, gx, gy);
                            if(i != gx && j != gy && bi != x && bj != y) {
                                return true;
                            }
                        }   
                    }
                }
            }
        }
    }

    return false;
}

void
blk_propagate(wfc_blocks_ptr blocks,
              uint32_t gx, uint32_t gy,
              uint64_t collapsed)
{
    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (int i = 0; i < blk_size; i++) {
        uint64_t idx = gx * blocks->grid_side * blk_size + gy * blk_size + i;
        if (bitfield_count(blocks->states[idx]) > 1) {
            blocks->states[idx] &= ~(collapsed);
        }
    }

    return;
}

void
grd_propagate_row(wfc_blocks_ptr blocks,
                  uint32_t gx, uint32_t gy, uint32_t x, uint32_t y,
                  uint64_t collapsed)
{
    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (int i = 0; i < blocks->grid_side; i++) {
        for (int j = 0; j < blocks->block_side; j++) {
            uint64_t idx = gx * blocks->grid_side * blk_size + i * blk_size + x * blocks->block_side + j;
            if (bitfield_count(blocks->states[idx]) > 1) {
                blocks->states[idx] &= ~(collapsed);
            }
        }
    }

    return;
}

void
grd_propagate_column(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy,
                     uint32_t x, uint32_t y, uint64_t collapsed)
{
    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (int i = 0; i < blocks->grid_side; i++) {
        for (int j = 0; j < blocks->block_side; j++) {
            uint64_t idx = i * blocks->grid_side * blk_size + gy * blk_size + j * blocks->block_side + y;
            if (bitfield_count(blocks->states[idx]) > 1) {
                blocks->states[idx] &= ~(collapsed);
            }
        }
    }

    return;
}

void printBinary2(uint64_t number) {
    // Determine the number of bits in u64_t
    int numBits = sizeof(uint64_t) + 1;

    // Loop through each bit in the number, starting from the most significant bit
    for (int i = numBits - 1; i >= 0; i--) {
        // Use a bitwise AND operation to check the value of the current bit
        if ((number & (1ULL << i)) != 0) {
            printf("1");
        } else {
            printf("0");
        }

        // Add space for better readability
        // if (i % 8 == 0) {
        //     printf(" ");
        // }
    }
    // printf("\n");
}

void grd_print(FILE *const file, const wfc_blocks_ptr block){
    FILE * fp = file;
    if(fp == NULL)
        fp = stdout;
    
    uint8_t gs = block->grid_side;
    uint8_t bs = block->block_side;

    for(uint32_t ii = 0; ii < gs; ii++){

        for(uint32_t i = 0; i < gs; i++){
            fprintf(fp, "+");
            for(uint32_t j = 0; j < bs; j++){
                fprintf(fp, "----------+");
            }
            fprintf(fp, "   ");
        }
        fprintf(fp, "\n");

        for(uint32_t jj = 0; jj < bs; jj++){

            for(uint32_t i = 0; i < gs; i++){
                fprintf(fp, "|");
                for(uint32_t j = 0; j < bs; j++){
                    const uint64_t collapsed = *blk_at(block, ii,i,jj,j);
                    printBinary2(collapsed);
                    fprintf(fp, " |", collapsed);
                    // fprintf(fp, "%3lu|", collapsed);
                }
                fprintf(fp, "   ");
            }
            fprintf(fp, "\n");

            for(int i = 0; i < gs; i++){
                fprintf(fp, "+");
                for(int j = 0; j < bs; j++){
                    fprintf(fp, "----------+");
                }
                fprintf(fp, "   ");
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }

}
