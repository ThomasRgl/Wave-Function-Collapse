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

    wfc_clone_DTD( blocks, init);

    uint32_t gs = blocks->grid_side;
    uint32_t bs = blocks->block_side;
    const uint64_t state_count = gs * gs * bs * bs;
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

    return ;

}

bool
solve_cuda(wfc_blocks_ptr d_blocks, wfc_blocks_ptr d_init, uint64_t seed)
{

    wfc_blocks buffer;
    cudaMemcpy(&buffer, d_init, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);

    // printf("solver: addr block : %p\n", blocks);

    dim3 dimGrid(1, 1, 1);
    dim3 dimBlock(buffer.block_side, buffer.block_side, 1);

    checkCudaErrors(cudaGetLastError());

    solve_cuda_device<<<dimGrid, dimBlock>>>(d_blocks, d_init, seed);

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    wfc_blocks * blocks =  wfc_clone_DTH( d_blocks);

    bool success = blocks->solved;
    if( success ){
        printf("Host : success with seed : %lu\n\n", blocks->seed);
        grd_print(NULL, blocks);
    }

    super_safe_free(blocks);
    return success;
}
