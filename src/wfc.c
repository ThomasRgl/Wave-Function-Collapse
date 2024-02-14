#include <stdint.h>
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
    
    return bitfield_count(state);
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
    uint64_t minima;
    uint8_t min_entropy = UINT8_MAX;
    vec2 the_location   = { 0 };
    uint64_t * blk_state = grd_at(blocks, gx, gy);

    uint32_t bs = blocks->block_side; 
    for (uint32_t y = 0; y < bs; y++) {
        for (uint32_t x = 0; x < bs; x++) {
            uint32_t idx = y * bs + x;
            uint64_t state = blk_state[ idx ] ;
            uint8_t  entropy = entropy_compute(state);

            // this condition accept 0
            // if entropy == 0, error will be raise later 
            if( entropy != 1  ){
                // smallest entropy
                if(min_entropy > entropy){
                    min_entropy = entropy;
                    minima = 0;
                    bitfield_set(minima, (uint8_t)idx);
                }
                //add to potential candidates
                else if (min_entropy == entropy) {
                    bitfield_set(minima, (uint8_t)idx);
                }
            }
        }
    }
    
    uint32_t candidate = bitfield_count(minima);
    entropy_location loc;
    loc.location = the_location;
    loc.entropy = min_entropy;
    return loc;
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
grd_check_error_in_column(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
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
grd_check_error_in_row(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
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
grd_check_error_in_block(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
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

bool
blk_propagate(wfc_blocks_ptr blocks,
              uint32_t gx, uint32_t gy,
              uint64_t collapsed, uint32_t * stack_cells, uint32_t * stack_size)
{
    bool changed = false;
    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (int i = 0; i < blk_size; i++) {
        uint64_t idx = gx * blocks->grid_side * blk_size + gy * blk_size + i;
        uint8_t entropy = entropy_compute(blocks->states[idx]);
        // if the cell is not collapsed yet
        if (entropy > 1) {
            blocks->states[idx] &= ~(collapsed);
            changed = true;
            // if the new entropy is 1 add the cell to the stack to propagate it later
            entropy = entropy_compute(blocks->states[idx]);
            if (entropy == 1) {
                stack_cells[*stack_size] = idx;
                (*stack_size)++;
            } else if (entropy == 0) {
                fprintf(stderr, "error in block propagation in block (%u, %u) at %u, %u\n", gx, gy, i % blocks->block_side, i / blocks->block_side);
                exit(EXIT_FAILURE);
            }
        }
    }

    return changed;
}

bool
grd_propagate_row(wfc_blocks_ptr blocks,
                  uint32_t gx, uint32_t gy, uint32_t x, uint32_t y,
                  uint64_t collapsed, uint32_t * stack_cells, uint32_t * stack_size)
{
    bool changed = false;
    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (int i = 0; i < blocks->grid_side; i++) {
        for (int j = 0; j < blocks->block_side; j++) {
            uint64_t idx = gx * blocks->grid_side * blk_size + i * blk_size + x * blocks->block_side + j;
            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                changed = true;
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    stack_cells[*stack_size] = idx;
                    (*stack_size)++;
                } else if (entropy == 0) {
                    fprintf(stderr, "error in row propagation in block (%u, %u) at %u, %u\n", gx, gy, x, y);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    return changed;
}

bool
grd_propagate_column(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy,
                     uint32_t x, uint32_t y, uint64_t collapsed, uint32_t * stack_cells, uint32_t * stack_size)
{
    bool changed = false;
    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (int i = 0; i < blocks->grid_side; i++) {
        for (int j = 0; j < blocks->block_side; j++) {
            uint64_t idx = i * blocks->grid_side * blk_size + gy * blk_size + j * blocks->block_side + y;
            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                changed = true;
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    stack_cells[*stack_size] = idx;
                    (*stack_size)++;
                } else if (entropy == 0) {
                    fprintf(stderr, "error in column propagation in block (%u, %u) at %u, %u\n", gx, gy, x, y);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    return changed;
}

bool
grd_propagate_all(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y, uint64_t collapsed)
{
    bool changed = false;

    uint32_t * stack_cells = malloc(sizeof(uint32_t) * 100);
    uint32_t stack_size = 0;

    // propagate the initial cell
    *blk_at(blocks, gx, gy, x, y) = collapsed;
    changed |= blk_propagate(blocks, gx, gy, collapsed, stack_cells, &stack_size);
    changed |= grd_propagate_column(blocks, gx, gy, x, y, collapsed, stack_cells, &stack_size);
    changed |= grd_propagate_row(blocks, gx, gy, x, y, collapsed, stack_cells, &stack_size);
    *blk_at(blocks, gx, gy, x, y) = collapsed;
    
    while (stack_size > 0) {

        // printf("stack (%d): ", stack_size);
        // for (int i = 0; i < stack_size; i++) {
        //     printf("%d ", stack_cells[i]);
        // }
        // printf("\n");

        // pop the cell from the stack
        uint64_t cell = stack_cells[0];
        for (uint32_t i = 0; i < stack_size; i++) {
            stack_cells[i] = stack_cells[i + 1];
        }
        stack_size--;

        // compute the coordinates of the new cell
        x = (cell / blocks->block_side) % blocks->block_side;
        y = cell % blocks->block_side;
        gx = cell / (blocks->grid_side * blocks->block_side * blocks->block_side);
        gy = (cell / (blocks->block_side * blocks->block_side)) % blocks->grid_side;

        // get the new collapsed state
        collapsed = *blk_at(blocks, gx, gy, x, y);

        // printf("collapsed (stack): %lu (", collapsed);
        // printBinary2(collapsed);
        // printf(") at : %lu, %lu, %lu, %lu\n", gx, gy, x, y);

        // propagate the new cell
        changed |= blk_propagate(blocks, gx, gy, collapsed, stack_cells, &stack_size);
        changed |= grd_propagate_column(blocks, gx, gy, x, y, collapsed, stack_cells, &stack_size);
        changed |= grd_propagate_row(blocks, gx, gy, x, y, collapsed, stack_cells, &stack_size);
        *blk_at(blocks, gx, gy, x, y) = collapsed;
    }

    free(stack_cells);

    return changed;
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
