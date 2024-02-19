#include "utils.cuh"
#include "types.cuh"
#include "bitfield.cuh"

#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strings.h>




wfc_blocks *
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



void
super_safe_free( wfc_blocks * blocks )
{
    free(blocks->stack_cells);
    free(blocks->row_masks);
    free(blocks->col_masks);
    free(blocks->blk_masks);
    free(blocks->states);
    free(blocks);

}

void
super_safe_Cudafree( wfc_blocks * d_blocks )
{
    wfc_blocks h_blocks;
    checkCudaErrors(cudaMemcpy(&h_blocks, d_blocks, sizeof(wfc_blocks), cudaMemcpyDeviceToHost));

    cudaFree(h_blocks.stack_cells);
    cudaFree(h_blocks.row_masks);
    cudaFree(h_blocks.col_masks);
    cudaFree(h_blocks.blk_masks);
    cudaFree(h_blocks.states);
    cudaFree(d_blocks);

}
// wfc_blocks *
// cloneToDevice( wfc_blocks * blocks, uint64_t seed)
// {
//     wfc_blocks * d_blocks;
//     cudaMalloc((void**)&d_blocks, sizeof(wfc_blocks));
//     // printf("clone: addr block : %p\n", d_blocks);

//     uint32_t gs = blocks->grid_side;
//     uint32_t bs = blocks->block_side;
//     const uint64_t state_count = gs * gs * bs * bs;

//     wfc_blocks * buffer = (wfc_blocks *) malloc(sizeof(wfc_blocks));
//     memcpy( buffer, blocks, sizeof(wfc_blocks));
//     buffer->seed = seed;
//     buffer->grid_side = gs;
//     buffer->block_side = bs;
//     buffer->solved = false;

//     checkCudaErrors(cudaMalloc((void**)&buffer->states    , state_count * sizeof(uint64_t) ));
//     checkCudaErrors(cudaMalloc((void**)&buffer->row_masks , gs * bs * sizeof(uint64_t) ));
//     checkCudaErrors(cudaMalloc((void**)&buffer->col_masks , gs * bs * sizeof(uint64_t) ));
//     checkCudaErrors(cudaMalloc((void**)&buffer->blk_masks , gs * gs * sizeof(uint64_t) ));
//     checkCudaErrors(cudaMalloc((void**)&buffer->stack_cells, (state_count-1) * sizeof(vec4) ));

//     checkCudaErrors(cudaMemcpy(buffer->states,    blocks->states,    state_count * sizeof(uint64_t), cudaMemcpyHostToDevice));
//     checkCudaErrors(cudaMemcpy(buffer->row_masks, blocks->row_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice));
//     checkCudaErrors(cudaMemcpy(buffer->col_masks, blocks->col_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice));
//     checkCudaErrors(cudaMemcpy(buffer->blk_masks, blocks->blk_masks, gs * gs * sizeof(uint64_t), cudaMemcpyHostToDevice));
    
//     checkCudaErrors(cudaMemcpy(d_blocks, buffer, sizeof(wfc_blocks), cudaMemcpyHostToDevice));

//     checkCudaErrors(cudaGetLastError());

//     free(buffer);
//     // printf("addr block : %p\n", d_blocks);

//     //
//     return d_blocks;    
// }

wfc_blocks * wfc_clone_HTD( wfc_blocks * src)
{
    wfc_blocks * d_dst;
    cudaMalloc((void**)&d_dst, sizeof(wfc_blocks));
    // printf("clone: addr block : %p\n", d_blocks);

    uint32_t gs = src->grid_side;
    uint32_t bs = src->block_side;
    const uint64_t state_count = gs * gs * bs * bs;

    wfc_blocks * buffer = (wfc_blocks *) malloc(sizeof(wfc_blocks));
    memcpy( buffer, src, sizeof(wfc_blocks));
    buffer->seed = 0;
    buffer->grid_side = (uint8_t)gs;
    buffer->block_side = (uint8_t)bs;
    buffer->solved = false;
    buffer->stack_size = 0;

    checkCudaErrors(cudaMalloc((void**)&buffer->states    , state_count * sizeof(uint64_t) ));
    checkCudaErrors(cudaMalloc((void**)&buffer->row_masks , gs * bs * sizeof(uint64_t) ));
    checkCudaErrors(cudaMalloc((void**)&buffer->col_masks , gs * bs * sizeof(uint64_t) ));
    checkCudaErrors(cudaMalloc((void**)&buffer->blk_masks , gs * gs * sizeof(uint64_t) ));
    checkCudaErrors(cudaMalloc((void**)&buffer->stack_cells, (state_count-1) * sizeof(vec4) ));

    checkCudaErrors(cudaMemcpy(buffer->states,    src->states,    state_count * sizeof(uint64_t), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(buffer->row_masks, src->row_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(buffer->col_masks, src->col_masks, gs * bs * sizeof(uint64_t), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(buffer->blk_masks, src->blk_masks, gs * gs * sizeof(uint64_t), cudaMemcpyHostToDevice));
    
    checkCudaErrors(cudaMemcpy(d_dst, buffer, sizeof(wfc_blocks), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaGetLastError());
    // printf("addr block : %p\n", d_blocks);

    free(buffer);
    //
    return d_dst;    
}


wfc_blocks * wfc_clone_DTH( wfc_blocks * d_src)
{

    wfc_blocks * buffer = (wfc_blocks *) malloc(sizeof(wfc_blocks));
    checkCudaErrors(cudaMemcpy(buffer, d_src, sizeof(wfc_blocks), cudaMemcpyDeviceToHost));

    uint32_t gs = buffer->grid_side;
    uint32_t bs = buffer->block_side;
    const uint64_t state_count = gs * gs * bs * bs;

    wfc_blocks * dst;
    dst = super_safe_malloc(gs, bs);

    checkCudaErrors(cudaMemcpy(dst->states,    buffer->states,    state_count * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(dst->row_masks, buffer->row_masks, gs * bs * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(dst->col_masks, buffer->col_masks, gs * bs * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(dst->blk_masks, buffer->blk_masks, gs * gs * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    dst->seed = buffer->seed;
    dst->solved = buffer->solved;
    dst->grid_side = (uint8_t)gs;
    dst->block_side = (uint8_t)bs;
    dst->stack_size = 0;
 
    checkCudaErrors(cudaGetLastError());

    free(buffer);
    //
    return dst;    
}

__device__
void wfc_clone_DTD( wfc_blocks * d_dst, wfc_blocks * d_src)
{

    const uint32_t gs = d_src->grid_side;
    const uint32_t bs = d_src->block_side;
    const uint64_t state_count = gs * gs * bs * bs;

    memcpy(d_dst->states,    d_src->states,    state_count * sizeof(uint64_t));
    memcpy(d_dst->row_masks, d_src->row_masks, gs * bs * sizeof(uint64_t));
    memcpy(d_dst->col_masks, d_src->col_masks, gs * bs * sizeof(uint64_t));
    memcpy(d_dst->blk_masks, d_src->blk_masks, gs * gs * sizeof(uint64_t));

    d_dst->seed = d_src->seed;
    d_dst->solved = d_src->solved;
    d_dst->grid_side = d_src->grid_side;
    d_dst->block_side = d_src->block_side;
    d_dst->stack_size = 0;
 

}
