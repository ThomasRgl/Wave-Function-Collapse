// #define _GNU_SOURCE

#include "bitfield.cuh"
#include "wfc.cuh"
#include <cstdio>

// #if !defined(WFC_CUDA)
// #error "WDC_CUDA should be defined..."
// #endif
__global__ void
solve_cuda_device(wfc_blocks_ptr blocks)
{
    printf("hello world from Device : %d \n", blocks->seed);
    return;
}

bool
solve_cuda(wfc_blocks_ptr blocks)
{
    solve_cuda_device<<<1,1>>>(NULL);
    return true;
}
