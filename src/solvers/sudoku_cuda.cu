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
solve_cuda_device(wfc_blocks_ptr blocks, wfc_blocks_ptr init, uint64_t seed)
{
    uint32_t gs = init->grid_side;
    uint32_t bs = init->block_side;
    const uint64_t state_count = gs * gs * bs * bs;

    blocks->block_side = bs; 
    blocks->grid_side = gs; 
    blocks->seed = seed; 
    memcpy(blocks->states,    init->states,    state_count * sizeof(uint64_t) );
    memcpy(blocks->row_masks, init->row_masks, gs * bs * sizeof(uint64_t) );
    memcpy(blocks->col_masks, init->col_masks, gs * bs * sizeof(uint64_t) );
    memcpy(blocks->blk_masks, init->blk_masks, gs * gs * sizeof(uint64_t) );

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


        bool error = grd_propagate_all(blocks, loc.gx,
                            loc.gy, loc.x, loc.y, collapsed_state);
        
        __syncthreads();

        if( error){
            break;
        }

        // if( error && threadIdx.x == 0 && threadIdx.y == 0){
        //     printf("error\n");
        // }

        iteration += 1;

        // if( threadIdx.x == 0 && threadIdx.y == 0){
        //     grd_print(NULL, blocks);
        // }
    }
    
    if(success && threadIdx.x == 0 && threadIdx.y == 0){
        // grd_print(NULL, blocks);
        blocks->solved = success; 
    }
    blocks->seed = 999;

    // getchar();
    return ;

}

bool
solve_cuda(wfc_blocks_ptr blocks, wfc_blocks_ptr d_init, uint64_t seed)
{

    // printf("solver: addr block : %p\n", blocks);

    dim3 dimGrid(1, 1, 1);
    dim3 dimBlock(3, 3, 1);

    checkCudaErrors(cudaGetLastError());

    solve_cuda_device<<<dimGrid, dimBlock>>>(blocks, d_init, seed);

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    // wfc_blocks_ptr b = (wfc_blocks*) malloc(sizeof(wfc_blocks));


    wfc_blocks_ptr b = super_safe_malloc( 3, 3);
    wfc_blocks_ptr buffer = super_safe_malloc( 3, 3);
    cudaMemcpy(buffer, blocks, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);
    cudaMemcpy(b->states, buffer->states, 3 * 3 * 3 * 3 * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    b->grid_side = 3;
    b->block_side = 3;
    bool succes = buffer->solved;
    if( succes ){
        printf("Host : success with seed : %lu\n\n", buffer->seed);
        grd_print(NULL, b);
    }
    else{
        printf("Host : failed\n");
    }

    // super_safe_free(b);
    // super_safe_free(buffer);
    return succes;
}
