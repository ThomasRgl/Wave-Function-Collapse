#include <stdint.h>
#define _GNU_SOURCE

#include "wfc.h"
#include "bitfield.h"
#include "utils.h" 
#include "md5.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <math.h>

// uint64_t 
// generate_random(
// uint32_t gx, uint32_t gy, uint32_t x, uint32_t y,
//                        uint64_t seed,
//                        uint64_t iteration
// ){}

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

    uint64_t state_count = bitfield_count(state);
    uint8_t random = (uint8_t)((uint64_t*) &digest)[0];
    uint8_t nth = random%state_count + 1;
    // printf(" nth state chosen : %u \n", (uint8_t) nth);
    // printf(" ancien state :  " );
    // printBinary2(state);
    // printf("  =  %lu", state);
    // printf("\n");
    // state = bitfield_only_nth_set(state, 2);
    // state = bitfield_only_nth_set(state, (uint8_t)random%state_count - 1);
    state = bitfield_only_nth_set(state, (uint8_t)nth - 1);
    // printf(" nouveau state : " );
    // printBinary2(state);
    // printf("  =  %lu", state);
    // printf("\n");

    return state;
}

uint8_t
entropy_compute(uint64_t state)
{
     
    return bitfield_count(state);
}

void
wfc_clone_into(wfc_blocks_ptr *const restrict ret_ptr, uint64_t seed, const wfc_blocks_ptr blocks)
{
    const uint64_t gs  = blocks->grid_side;
    const uint64_t bs = blocks->block_side;
    wfc_blocks_ptr ret        = *ret_ptr;

    // TODO
    if (NULL == ret) {
        if (NULL == (ret = super_safe_malloc(gs, bs)) ) {
            fprintf(stderr, "failed to clone blocks structure\n");
            exit(EXIT_FAILURE);
        }
    } else if ( gs != ret->grid_side || bs != ret->block_side) {
        fprintf(stderr, "size mismatch!\n");
        exit(EXIT_FAILURE);
    }
    
    
    memcpy(ret->states, blocks->states, gs * gs * bs * bs * sizeof(uint64_t) );
    memcpy(ret->row_masks, blocks->row_masks, gs * bs * sizeof(uint64_t) );
    memcpy(ret->col_masks, blocks->col_masks, gs * bs * sizeof(uint64_t) );
    memcpy(ret->blk_masks, blocks->blk_masks, gs * gs * sizeof(uint64_t) );
    ret->seed = seed;
    ret->block_side = (uint8_t)bs;
    ret->grid_side = (uint8_t)gs;

    *ret_ptr       = ret;
}

// uint64_t *
// choose_and_collapse (const wfc_blocks_ptr blocks, uint64_t ite){
//     uint8_t min_entropy = UINT8_MAX;
//     uint32_t loc_x = 0 ; 
//     uint32_t loc_y = 0;
//     uint32_t loc_gy = 0;
//     uint32_t loc_gx = 0;
//
//     for( uint32_t gy = 0; gy < blocks->grid_side; gy++ ){
//         for( uint32_t gx = 0; gx < blocks->grid_side; gx++ ){
//             entropy_location loc = blk_min_entropy(blocks, gx, gy);
//             if( loc.entropy < min_entropy ){
//                 min_entropy = loc.entropy;
//                 loc_x = loc.location.x;
//                 loc_y = loc.location.y;
//                 loc_gy = gy;
//                 loc_gx = gx;
//             }
//         }
//     }
//     printf(" choose loc : [%d, %d] [%d, %d] = %d\n", loc_gy, loc_gx, loc_y, loc_x, min_entropy );
//     
//     uint64_t * state = blk_at(blocks, loc_gx, loc_gy, loc_x, loc_y);
//     uint64_t collapsed_state = entropy_collapse_state(
//         *state, loc_gx, loc_gy, loc_x, loc_y, blocks->seed, ite);
//
//     *state = collapsed_state;
//
//     blk_propagate(blocks, loc_gx, loc_gy, collapsed_state);
//     grd_propagate_column(blocks, loc_gx, loc_gy, loc_x, loc_y, collapsed_state);
//     grd_propagate_row(blocks, loc_gx, loc_gy, loc_x, loc_y, collapsed_state);
//  
//
//     return state;
//     
// }

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
            uint64_t idx = idx_at(blocks, 0, 0, x, y);
            uint64_t state = blk_state[ idx ] ;
            uint8_t  entropy = entropy_compute(state);
            // printf("idx %d %d = %d == %d\n",y , x, idx, entropy);

            // this condition accept 0
            // if entropy == 0, error will be raise later 
            if( entropy != 1  ){
                // smallest entropy
                if(min_entropy > entropy){
                    min_entropy = entropy;
                    the_location.y = y;
                    the_location.x = x;

                    // minima = 0; // minima store the location of each minimal and equal state,
                    // when we find a smaller entropy we reset that bloc
                    // bitfield_set(minima, (uint8_t)idx);
                }
                //add to potential candidates
                // else if (min_entropy == entropy) {
                //     bitfield_set(minima, (uint8_t)idx);
                // }
            }
        }
    }
    
    // uint32_t candidate = bitfield_count(minima);
    entropy_location loc;
    loc.location = the_location;
    loc.entropy = min_entropy;

    // printf("min entropy blk (%u; %u) : [%u; %u] = %u\n",
    //                gy, gx, loc.location.y, loc.location.x, loc.entropy );

    return loc;
}

// static inline uint64_t
// blk_filter_mask_for_row(wfc_blocks_ptr blocks,
//                            uint32_t gx, uint32_t gy, uint32_t y)
// {
//     uint64_t mask = 0;
//     uint64_t state = 0;
//     for (uint32_t i = 0; i < blocks->block_side; i++) {
//         state = *blk_at(blocks, gx, gy, i, y);
//         if(bitfield_count(state) == 1){
//             mask |= state;
//         }
//     }
//
//     return mask;
// }
//
// static inline uint64_t
// blk_filter_mask_for_column(wfc_blocks_ptr blocks,
//                             uint32_t gx, uint32_t gy, uint32_t x)
// {
//     uint64_t mask = 0;
//     uint64_t state = 0;
//     for (uint32_t i = 0; i < blocks->block_side; i++) {
//         state = *blk_at(blocks, gx, gy, x, i);
//         if(bitfield_count(state) == 1){
//             mask |= state;
//         }
//     }
//
//     return mask;
// }
//
// static inline uint64_t
// grd_filter_mask_for_row(wfc_blocks_ptr blocks,
//                            uint32_t gy, uint32_t y)
// {
//     uint64_t mask = 0;
//     for (uint32_t i = 0; i < blocks->grid_side; i++) {
//         mask |= blk_filter_mask_for_row(blocks, i, gy, y);
//     }
//     
//     return mask;
// }
//
// static inline uint64_t
// grd_filter_mask_for_column(wfc_blocks_ptr blocks,
//                         uint32_t gx, uint32_t x)
// {
//     uint64_t mask = 0;
//     for (uint32_t i = 0; i < blocks->grid_side; i++) {
//         mask |= blk_filter_mask_for_column(blocks, gx, i, x);
//     }
//
//     return mask;
// }
//
// static inline uint64_t
// blk_filter_mask_for_block(wfc_blocks_ptr blocks,
//                           uint32_t gy, uint32_t gx)
// {
//     uint64_t mask = 0;
//     uint64_t state = 0;
//     for (uint32_t i = 0; i < blocks->block_side; i++) {
//         for (uint32_t j = 0; j < blocks->block_side; j++) {
//             state = *blk_at(blocks, gx, gy, i, j);
//             if(bitfield_count(state) == 1){
//                 mask |= state;
//             }
//         }
//     }
//
//     return mask;
// }
//
// bool
// grd_check_error_in_column(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
// {
//     uint64_t mask = grd_filter_mask_for_column(blocks, gx, x);
//     uint64_t state = 0;
//     for (uint32_t i = 0; i < blocks->grid_side; i++) {
//         for (uint32_t j = 0; j < blocks->block_side; j++) {
//             state = *blk_at(blocks, gx, i, x, j);
//             if (bitfield_count(state) == 1) {
//                 if (bitfield_count(state|mask) != bitfield_count(mask)+1) {
//                     // printf("Col : i: %u, j: %u, x: %u, y: %u, gx: %u, gy: %u\n", i, j, x, y, gx, gy);
//                     if(i != gy && j != y) {
//                         return true;
//                     }
//                 }   
//             }
//         }
//     }
//   
//     return false;
// }
//
// bool
// grd_check_error_in_row(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
// {
//     uint64_t mask = grd_filter_mask_for_row(blocks, gy, y);
//     uint64_t state = 0;
//     for (uint32_t i = 0; i < blocks->grid_side; i++) {
//         for (uint32_t j = 0; j < blocks->block_side; j++) {
//             state = *blk_at(blocks, i, gy, j, y);
//             if (bitfield_count(state) == 1) {
//                 if (bitfield_count(state|mask) != bitfield_count(mask)+1) {
//                     // printf("Row : i: %u, j: %u, x: %u, y: %u, gx: %u, gy: %u\n", i, j, x, y, gx, gy);
//                     if(i != gx && j != x) {
//                         return true;
//                     }
//                 }   
//             }
//         }
//     }
//
//     return false;
// }
//
// bool
// grd_check_error_in_block(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
// {
//     uint64_t mask = blk_filter_mask_for_block(blocks, gx, gy);
//     uint64_t state = 0;
//     for (uint32_t i = 0; i < blocks->grid_side; i++) {
//         for (uint32_t j = 0; j < blocks->grid_side; j++) {
//             for (uint32_t bi = 0; bi < blocks->block_side; bi++) {
//                 for (uint32_t bj = 0; bj < blocks->block_side; bj++) {
//                     state = *blk_at(blocks, i, j, bi, bj);
//                     if (bitfield_count(state) == 1) {
//                         if (bitfield_count(state|mask) != bitfield_count(mask)+1) {
//                             // printf("Block : i: %u, j: %u, bi: %u, bj: %u, x: %u, y: %u, gx: %u, gy: %u\n", i, j, bi, bj, x, y, gx, gy);
//                             if(i != gx && j != gy && bi != x && bj != y) {
//                                 return true;
//                             }
//                         }   
//                     }
//                 }
//             }
//         }
//     }
//
//     return false;
// }

bool
blk_propagate(wfc_blocks_ptr blocks,
              uint32_t gx, uint32_t gy,
              uint64_t collapsed)
{

    if ( (blocks->blk_masks[gy * blocks->block_side + gx] & collapsed) == 0) {
        // fprintf(stderr, "error in (mask) block propagation in block (%u, %u)\n", gy, gx);
        return true;
    }

    blocks->blk_masks[gy * blocks->block_side + gx] &= ~collapsed;

    for (uint32_t y = 0; y < blocks->block_side; y++) {
        for (uint32_t x = 0; x < blocks->block_side; x++) {

            uint64_t idx = idx_at(blocks, gx, gy, x, y);
            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    vec4 coord = {gx, gy, x, y};
                    blocks->stack_cells[blocks->stack_size] = coord;
                    blocks->stack_size++;
                }
            }
        }
    }

    return true;
}

bool
grd_propagate_column(wfc_blocks_ptr blocks,
                  uint32_t gx, uint32_t __, uint32_t x, uint32_t _,
                  uint64_t collapsed)
{

    if ( (blocks->col_masks[gx * blocks->block_side + x] & collapsed) == 0) {
        // fprintf(stderr, "error in (mask) column propagation in column block %u in column %u\n", gx, x);
        return false;
    }

    blocks->col_masks[gx * blocks->block_side + x] &= ~collapsed;

    for (uint32_t gy = 0; gy < blocks->grid_side; gy++) {
        for (uint32_t y = 0; y < blocks->block_side; y++) {
            uint64_t idx = idx_at(blocks, gx, gy, x, y);

            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    vec4 coord = {gx, gy, x, y};
                    blocks->stack_cells[blocks->stack_size] = coord;
                    blocks->stack_size++;
                }
            }
        }
    }

    return true;
}


bool
grd_propagate_row(wfc_blocks_ptr blocks, uint32_t __, uint32_t gy,
                     uint32_t _, uint32_t y, uint64_t collapsed)
{

    if ( (blocks->row_masks[gy * blocks->block_side + y] & collapsed) == 0) {
        // fprintf(stderr, "error in (mask) row propagation in row block %u in row %u\n", gy, y);
        return false;
    }

    blocks->row_masks[gy * blocks->block_side + y] &= ~collapsed;

    uint32_t blk_size = blocks->block_side*blocks->block_side;
    for (uint32_t gx = 0; gx < blocks->grid_side; gx++) {
        for (uint32_t x = 0; x < blocks->block_side; x++) {
            // uint64_t idx = i * blocks->grid_side * blk_size + gy * blk_size + j * blocks->block_side + y;
            uint64_t idx = idx_at(blocks, gx, gy, x, y);
            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    vec4 coord = {gx, gy, x, y};
                    blocks->stack_cells[blocks->stack_size] = coord;
                    blocks->stack_size++;
                }
            }
        }
    }

    return true;
}

bool
grd_propagate_all(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y, uint64_t collapsed)
{
    bool no_error = true;

    // propagate the initial cell
    no_error &= blk_propagate(blocks, gx, gy, collapsed);
    no_error &= grd_propagate_column(blocks, gx, gy, x, y, collapsed);
    no_error &= grd_propagate_row(blocks, gx, gy, x, y, collapsed);
    *blk_at(blocks, gx, gy, x, y) = collapsed;
    
    while (no_error && blocks->stack_size > 0) {

        // printf("stack (%d): ", blocks->stack_size);
        // for (int i = 0; i < blocks->stack_size; i++) {
        //     printf("%d ", blocks->stack_cells[i]);
        // }
        // printf("\n");

        // pop the cell from the stack
        vec4 coord = blocks->stack_cells[0];
        x = coord.x;
        y = coord.y;
        gx = coord.gx;
        gy = coord.gy;
        for (uint32_t i = 0; i < blocks->stack_size; i++) {
            blocks->stack_cells[i] = blocks->stack_cells[i + 1];
        }
        blocks->stack_size--;

        // get the new collapsed state
        collapsed = *blk_at(blocks, gx, gy, x, y);

        // printf("collapsed (stack): %lu (", collapsed);
        // printBinary2(collapsed);
        // printf(") at : %lu, %lu, %lu, %lu\n", gy, gx, y, x);

        // propagate the new cell
        no_error &= blk_propagate(blocks, gx, gy, collapsed);
        no_error &= grd_propagate_column(blocks, gx, gy, x, y, collapsed);
        no_error &= grd_propagate_row(blocks, gx, gy, x, y, collapsed);
        *blk_at(blocks, gx, gy, x, y) = collapsed;
    }

    return no_error;
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
                    const uint64_t collapsed = *blk_at(block, i,ii,j,jj);
                    if( bitfield_count(collapsed) == 1)
                        fprintf(fp, "   %lu     ", (uint64_t)log2((double)collapsed)+1 );
                    else
                        printBinary2(collapsed);
                    // fprintf(fp, "   %lu     ", idx_at(block, i,ii,j,jj) );
                    fprintf(fp, " |");
                    // fprintf(fp, "%3lu|");
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
