#pragma once

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


#define getLastCudaError(msg) __getLastCudaError(msg, __FILE__, __LINE__)
#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }

static inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=false)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

static inline void __getLastCudaError(const char *errorMessage, const char *file,
                               const int line)
{
  cudaError_t err = cudaGetLastError();

  if (cudaSuccess != err)
  {
    fprintf(stderr,
            "%s(%i) : getLastCudaError() CUDA error :"
            " %s : (%d) %s.\n",
            file, line, errorMessage, (int)(err),
            cudaGetErrorString(err));
    exit(EXIT_FAILURE);
  }
}

wfc_blocks * super_safe_malloc( uint32_t gs, uint32_t bs);

void super_safe_free( wfc_blocks * blocks );
void super_safe_Cudafree( wfc_blocks * d_blocks );
wfc_blocks * cloneToDevice( wfc_blocks * blocks, uint64_t seed);


wfc_blocks * wfc_clone_HTD( wfc_blocks * src);
wfc_blocks * wfc_clone_DTH( wfc_blocks * d_src);

void wfc_clone_DTD( wfc_blocks * d_dst, wfc_blocks * d_src, uint32_t gs, uint32_t bs);
wfc_blocks * wfc_clone_HTH( wfc_blocks * dst, wfc_blocks * src);


