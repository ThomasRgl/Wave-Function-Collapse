// #define _GNU_SOURCE

#include "bitfield.cuh"
#include "wfc.cuh"
#include "utils.cuh"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>


// #if !defined(WFC_CUDA)
// #error "WDC_CUDA should be defined..."
// #endif
__global__ void
solve_cuda_device(wfc_blocks_ptr ret_blocks, wfc_blocks_ptr init, uint64_t seed)
{

    // extern __shared__ uint64_t big_array[];
    __shared__ uint64_t big_array[2999];

    __shared__ wfc_blocks b;
    wfc_blocks_ptr blocks = &b;

    if( threadIdx.x == 0 && threadIdx.y == 0){
        const uint32_t gs = init->grid_side;
        const uint32_t bs = init->block_side;
        const uint64_t state_count = gs * gs * bs * bs;

        blocks->states      = (uint64_t*) &big_array[0];
        blocks->row_masks   = blocks->states + state_count;
        blocks->col_masks   = blocks->row_masks + gs*bs;
        blocks->blk_masks   = blocks->col_masks + gs*bs;
        blocks->stack_cells = (vec4 *)(blocks->blk_masks + gs * gs);
    
        wfc_clone_DTD(blocks, init);
    }
    __syncthreads();




    // uint32_t gs = blocks->grid_side;
    // uint32_t bs = blocks->block_side;
    // const uint64_t state_count = gs * gs * bs * bs;
    blocks->seed = seed;

    uint64_t iteration  = 0;

    bool success = false;

    jusqua_la_retraite {

        // if( threadIdx.x == 0 && threadIdx.y == 0){
        //     printf("\n======== iteration : %lu ========\n\n", iteration);
        // }

        vec4 loc = grd_min_entropy(blocks);
        
        if( loc.x == UINT8_MAX ){
            success = true;
            // printf("success\n");
            break;
        }

        uint64_t * state = blk_at(blocks, loc.gx, loc.gy,
                                      loc.x, loc.y);
        
        if( state == 0){
            break;
        }

        // if( state == 0 && threadIdx.x == 0 && threadIdx.y == 0) {
        //     printf("state = 0\n");
        // }

        __shared__ uint64_t collapsed_state;

        if( threadIdx.x == 0 && threadIdx.y == 0){
            collapsed_state = entropy_collapse_state(
            *state, loc.gx, loc.gy, loc.x, loc.y,
            blocks->seed, iteration);
            *state = collapsed_state;

            // printf("collapsed (entropy): %lu (",(uint64_t)log2((double)collapsed_state)+1);
            // printBinary2(collapsed_state);
            // printf(") at : [%u, %u] [%u, %u]\n", loc.gy, loc.gx, loc.y, loc.x);
        }
        __syncthreads();

        bool error = grd_propagate_all(blocks, loc.gx,
                            loc.gy, loc.x, loc.y, collapsed_state);
        
        __syncthreads();

        

        // if( error && threadIdx.x == 0 && threadIdx.y == 0){
        //     printf("error\n");
        // }
    
        if( error){
            break;
        }

        iteration += 1;

        // if( threadIdx.x == 0 && threadIdx.y == 1){
        //     grd_print(NULL, blocks);
        // }
    }
        
    if(success && threadIdx.x == 0 && threadIdx.y == 0){
        // grd_print(NULL, blocks);
        blocks->solved = success; 
        wfc_clone_DTD(ret_blocks, blocks);
    }

    return ;

}

bool
solve_cuda(wfc_blocks_ptr d_blocks, wfc_blocks_ptr d_init, uint64_t seed)
{

    // wfc_blocks buffer;
    // cudaMemcpy(&buffer, d_init, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);

    // printf("solver: addr block : %p\n", blocks);

    dim3 dimGrid(1, 1, 1);
    // dim3 dimBlock(buffer.block_side, buffer.block_side, 1);
    dim3 dimBlock(6, 6, 1);

    checkCudaErrors(cudaGetLastError());

    const uint32_t gs = 6;
    const uint32_t bs = 6;
    const uint64_t state_count = gs * gs * bs * bs;
    
    uint32_t sm = 
        state_count + // states
        gs * bs +// row
        gs * bs +// col
        gs * gs +//blk
        state_count;//stack


    
    solve_cuda_device<<<dimGrid, dimBlock, 0 >>>(d_blocks, d_init, seed);

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    wfc_blocks * blocks =  wfc_clone_DTH( d_blocks); // can be optimized

    bool success = blocks->solved;
    if( success ){
        printf("Host : success with seed : %lu\n\n", blocks->seed);
        grd_print(NULL, blocks);
        verify_block(blocks);
    }

    super_safe_free(blocks);
    return success;
}
