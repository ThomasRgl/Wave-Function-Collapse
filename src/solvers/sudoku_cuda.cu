// #define _GNU_SOURCE

#include "bitfield.cuh"
#include "wfc.cuh"
#include "utils.cuh"

#include <cstdio>
#include <cstdlib>


// #if !defined(WFC_CUDA)
// #error "WDC_CUDA should be defined..."
// #endif
__global__ void
solve_cuda_device(wfc_blocks_ptr blocks)
{
    uint64_t iteration  = 0;
    const uint64_t seed = blocks->seed;
    // struct {
    //     uint32_t gy, x, y, _1;
    //     uint64_t state;
    // } row_changes[blocks->grid_side];

    // grd_print(NULL, blocks);
    // getchar();

    bool success = false;

    jusqua_la_retraite {

        // if( threadIdx.x == 0 && threadIdx.y == 0){
        //     printf("\n======== iteration : %lu ========\n\n", iteration);
        // }

        vec4 loc;
        
        loc = grd_min_entropy(blocks);
      
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
        grd_print(NULL, blocks);
        blocks->solved = success; 
    }
    // getchar();
    return ;
}

bool
solve_cuda(wfc_blocks_ptr blocks)
{

    // printf("solver: addr block : %p\n", blocks);

    dim3 dimGrid(1, 1, 1);
    dim3 dimBlock(3, 3, 1);

    checkCudaErrors(cudaGetLastError());
    solve_cuda_device<<<dimGrid, dimBlock>>>(blocks);

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    // wfc_blocks_ptr b = (wfc_blocks*) malloc(sizeof(wfc_blocks));

    wfc_blocks_ptr b = super_safe_malloc( 3, 3);
    wfc_blocks_ptr buffer = super_safe_malloc( 3, 3);
    cudaMemcpy(buffer, blocks, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);
    cudaMemcpy(b->states, buffer->states, 3 * 3 * 3 * 3 * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    printf("hello world from Host, success ? : %u \n", b->solved);
    b->grid_side = 3;
    b->block_side = 3;
    grd_print(NULL, b);
    bool succes = b->solved;
    // super_safe_free(b);
    // super_safe_free(buffer);
    return succes;
}
