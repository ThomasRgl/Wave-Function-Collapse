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

        // grd_print(NULL, blocks);

        vec4 loc = grd_min_entropy(blocks);
        // printf(" choose loc :   [%d, %d] : [%d, %d] = %lu\n", 
        // loc.gy, loc.gx, loc.y, loc.x, *state );
      
        if( loc.x == UINT8_MAX ){
            success = true;
            // printf("success\n");
            break;
        }

        uint64_t * state = blk_at(blocks, loc.gx, loc.gy,
                                      loc.x, loc.y);
 

        if( state == 0){
            // printf("state = 0\n");
            break;
        }

        uint64_t collapsed_state = entropy_collapse_state(
            *state, loc.gx, loc.gy, loc.x, loc.y,
            blocks->seed, iteration);
        *state = collapsed_state;

        bool error = grd_propagate_all(blocks, loc.gx,
                            loc.gy, loc.x, loc.y, collapsed_state);
        

        if( error ){
            // printf("error\n");
            break;
        }

        iteration += 1;
    }
    
    if(success){
	grd_print(NULL, blocks);
    }
    // getchar();
    blocks->solved = success; 
    return ;
}

bool
solve_cuda(wfc_blocks_ptr blocks)
{

    printf("solver: addr block : %p\n", blocks);

    checkCudaErrors(cudaGetLastError());
    solve_cuda_device<<<1,1>>>(blocks);

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    wfc_blocks_ptr b = (wfc_blocks*) malloc(sizeof(wfc_blocks));
    cudaMemcpy(b, blocks, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);
    printf("hello world from Host, seccess ? : %u \n", b->solved);
    // grd_print(NULL, b);
    bool succes = b->solved;
    free(b);
    return succes;
}
