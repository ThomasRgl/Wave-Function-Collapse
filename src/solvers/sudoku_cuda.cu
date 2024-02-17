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
    printf("hello world from Device : %lu \n", blocks->seed);
    grd_print(NULL, blocks);
    return;
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
    b->seed = 10;
    cudaMemcpy(b, blocks, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);
    printf("hello world from Host : %lu \n", b->seed);
    // grd_print(NULL, b);
    return true;
}
