#pragma once

#include "types.cuh"
#include "bitfield.cuh"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strings.h>

static inline wfc_blocks *
super_safe_malloc( uint32_t gs, uint32_t bs)
{

    const uint64_t state_count = gs * gs * bs * bs;
    wfc_blocks *ret = (wfc_blocks *)malloc( sizeof(wfc_blocks) );
    ret->states    = (uint64_t*) malloc( state_count * sizeof(uint64_t) );
    ret->row_masks = (uint64_t*) malloc( gs * bs * sizeof(uint64_t) );
    ret->col_masks = (uint64_t*) malloc( gs * bs * sizeof(uint64_t) );
    ret->blk_masks = (uint64_t*) malloc( gs * gs * sizeof(uint64_t) );
    ret->stack_cells = (vec4*) malloc( (state_count-1) * sizeof(vec4) );

    ret->stack_size = 0;

    uint64_t mask = 0;
    for (uint8_t i = 0; i < bs * bs; i += 1) {
        mask = bitfield_set(mask, i);
    }
    
    for (uint64_t i = 0; i < state_count; i += 1) 
        ret->states[i] = mask;

    for (uint64_t i = 0; i < gs * bs; i += 1) {
        ret->row_masks[i] = mask;
        ret->col_masks[i] = mask;
    }

    for (uint64_t i = 0; i < gs * gs; i += 1) 
        ret->blk_masks[i] = mask;

    return ret; 
}

// static inline wfc_blocks *
// super_safe_CudaMalloc( uint32_t gs, uint32_t bs)
// {
//
//     wfc_blocks *ret_d; 
//     const uint64_t state_count = gs * gs * bs * bs;
//     cudaMalloc((void**)&ret_d, sizeof(wfc_blocks));
//     cudaMalloc((void**)&ret_d->states, sizeof(wfc_blocks));
//     cudaMalloc((void**)ret_d->states    , state_count * sizeof(uint64_t) );
//     cudaMalloc((void**)ret_d->row_masks , gs * bs * sizeof(uint64_t) );
//     cudaMalloc((void**)ret_d->col_masks , gs * bs * sizeof(uint64_t) );
//     cudaMalloc((void**)ret_d->blk_masks , gs * gs * sizeof(uint64_t) );
//     cudaMalloc((void**)ret_d->stack_cells, (state_count-1) * sizeof(vec4) );
//
//     // ret->stack_cells = (vec4*) malloc( (state_count-1) * sizeof(vec4) );
//     // cudaMemcpy(ret_d->stack_size, , size_t count, enum cudaMemcpyKind kind)
//     // ret->stack_size = 0; // later
//
//     // uint64_t mask = 0;
//     // for (uint8_t i = 0; i < bs * bs; i += 1) {
//     //     mask = bitfield_set(mask, i);
//     // }
//     // 
//     // for (uint64_t i = 0; i < state_count; i += 1) 
//     //     ret->states[i] = mask;
//     //
//     // for (uint64_t i = 0; i < gs * bs; i += 1) {
//     //     ret->row_masks[i] = mask;
//     //     ret->col_masks[i] = mask;
//     // }
//     //
//     // for (uint64_t i = 0; i < gs * gs; i += 1) 
//     //     ret->blk_masks[i] = mask;
//
//     return ret_d; 
// }




static inline void
super_safe_free( wfc_blocks * blocks )
{
    free(blocks->stack_cells);
    free(blocks->row_masks);
    free(blocks->col_masks);
    free(blocks->blk_masks);
    free(blocks->states);
    free(blocks);

}

static inline void
super_safe_Cudafree( wfc_blocks * d_blocks )
{
    free(d_blocks->stack_cells);
    free(d_blocks->row_masks);
    free(d_blocks->col_masks);
    free(d_blocks->blk_masks);
    free(d_blocks->states);
    free(d_blocks);

}
static inline wfc_blocks *
cudaCloneToDevice( wfc_blocks * blocks)
{
    wfc_blocks * d_blocks;
    cudaMalloc((void**)&d_blocks, sizeof(wfc_blocks));

    uint32_t gs = blocks->grid_side;
    uint32_t bs = blocks->block_side;
    const uint64_t state_count = gs * gs * bs * bs;

    wfc_blocks * buffer = (wfc_blocks *) malloc(sizeof(wfc_blocks));
    memcpy( buffer, blocks, sizeof(wfc_blocks));
    cudaMalloc((void**)&buffer->states    , state_count * sizeof(uint64_t) );
    cudaMalloc((void**)&buffer->row_masks , gs * bs * sizeof(uint64_t) );
    cudaMalloc((void**)&buffer->col_masks , gs * bs * sizeof(uint64_t) );
    cudaMalloc((void**)&buffer->blk_masks , gs * gs * sizeof(uint64_t) );
    cudaMalloc((void**)&buffer->stack_cells, (state_count-1) * sizeof(vec4) );



    cudaMemcpy(d_blocks, buffer, sizeof(wfc_blocks), cudaMemcpyHostToDevice);

    //    cudaMemcpy(d_blocks->states,    blocks->states,    state_count * sizeof(uint64_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_blocks->row_masks, blocks->row_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_blocks->col_masks, blocks->col_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_blocks->blk_masks, blocks->blk_masks, gs * gs * sizeof(uint64_t), cudaMemcpyHostToDevice);

    return d_blocks;    
}

// static inline wfc_blocks *
// cudaCloneToDevice( wfc_blocks * blocks)
// {
//     wfc_blocks * d_blocks;
//
//     uint32_t gs = blocks->grid_side;
//     uint32_t bs = blocks->block_side;
//     const uint64_t state_count = gs * gs * bs * bs;
//
//     cudaMalloc((void**)&d_blocks, sizeof(wfc_blocks));
//     cudaMemcpy(d_blocks, blocks, sizeof(wfc_blocks), cudaMemcpyHostToDevice);
//
//     cudaMalloc((void**)&d_blocks->states    , state_count * sizeof(uint64_t) );
//     cudaMalloc((void**)&d_blocks->row_masks , gs * bs * sizeof(uint64_t) );
//     cudaMalloc((void**)&d_blocks->col_masks , gs * bs * sizeof(uint64_t) );
//     cudaMalloc((void**)&d_blocks->blk_masks , gs * gs * sizeof(uint64_t) );
//     cudaMalloc((void**)&d_blocks->stack_cells, (state_count-1) * sizeof(vec4) );
//
//     cudaMemcpy(d_blocks->states,    blocks->states,    state_count * sizeof(uint64_t), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_blocks->row_masks, blocks->row_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_blocks->col_masks, blocks->col_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_blocks->blk_masks, blocks->blk_masks, gs * gs * sizeof(uint64_t), cudaMemcpyHostToDevice);
//
//     return d_blocks;    
// }
