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
solve_cuda_device(wfc_blocks_ptr ret_blocks, wfc_blocks_ptr init, uint64_t * seed_list, uint32_t * solved)
{

    uint64_t seed = seed_list[blockIdx.x];

    // extern __shared__ uint64_t big_array[];
    __shared__ uint64_t big_array[4999];

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


    // wfc_blocks_ptr blocks = ret_blocks;
    // wfc_clone_DTD(blocks, init);



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
        //     // grd_print(NULL, blocks);
        // }
    
        if( error){
            break;
        }

        iteration += 1;

        // if( threadIdx.x == 0 && threadIdx.y == 0){
        //     grd_print(NULL, blocks);
        // }
    }
        
    if(success && threadIdx.x == 0 && threadIdx.y == 0){
        blocks->solved = success; 
    }

    // if(threadIdx.x == 0 && threadIdx.y == 0){
    //     grd_print(NULL, blocks);
    // }

    ////////////////////////////////////
    uint32_t first = 1; 
    if(success && threadIdx.x == 0 && threadIdx.y == 0){
        first = atomicAdd(solved, 1);
        // printf("B%u A trouvé !  val : %u  ;  SEED = %lu\n", blockIdx.x, first, seed);

        if( first == 0 ){
            grd_print(NULL, blocks);
            wfc_clone_DTD(ret_blocks, blocks);
        }
    }
    // else if( threadIdx.x == 0 && threadIdx.y == 0){
    //     printf("B%u n'a PAS trouvé !  SEED : %lu \n", blockIdx.x, seed);
    // }

    return ;

}

bool
solve_cuda(wfc_blocks_ptr d_blocks, wfc_blocks_ptr d_init, uint64_t * seed_list, uint32_t gs, uint32_t bs, uint64_t p)
{

    // wfc_blocks buffer;
    // cudaMemcpy(&buffer, d_init, sizeof(wfc_blocks), cudaMemcpyDeviceToHost);

    // printf("solver: addr block : %p\n", blocks);

    dim3 dimGrid(p, 1, 1);
    // dim3 dimBlock(buffer.block_side, buffer.block_side, 1);
    dim3 dimBlock(gs, gs, 1);

    checkCudaErrors(cudaGetLastError());

    const uint64_t state_count = gs * gs * bs * bs;
    
    uint64_t sm = 
        state_count + // states
        gs * bs +// row
        gs * bs +// col
        gs * gs +//blk
        state_count;//stack


    uint64_t * d_seed_list;
    checkCudaErrors(cudaMalloc((void**)&d_seed_list, p * sizeof(uint64_t) ));
    checkCudaErrors(cudaMemcpy(d_seed_list, seed_list, p * sizeof(uint64_t), cudaMemcpyHostToDevice));

    uint32_t solved = 0;
    uint32_t * d_solved;
    checkCudaErrors(cudaMalloc((void**)&d_solved, sizeof(uint32_t) ));
    checkCudaErrors(cudaMemcpy(d_solved, &solved, sizeof(uint32_t), cudaMemcpyHostToDevice));



    checkCudaErrors(cudaGetLastError());

    // uint64_t is_solved = 0;
    // uint64_t * d_solved = 0;
    // checkCudaErrors(cudaMalloc((void**)&d_solved, sizeof(uint64_t) ));

    // checkCudaErrors(cudaMalloc((void**)&buffer->states    , state_count * sizeof(uint64_t) ));
    
    solve_cuda_device<<<dimGrid, dimBlock, 0 >>>(d_blocks, d_init, d_seed_list, d_solved);

    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());


    checkCudaErrors(cudaMemcpy(&solved, d_solved, sizeof(uint32_t), cudaMemcpyDeviceToHost));

   

    if( solved >= 1){
        wfc_blocks * blocks =  wfc_clone_DTH( d_blocks); // can be optimized
        printf("Host : success with seed : %lu\n\n", blocks->seed);
        // grd_print(NULL, blocks);
        verify_block(blocks);
        super_safe_free(blocks);
    }

    cudaFree(d_seed_list);
    cudaFree(d_solved);



    return ( solved >= 1 ? true : false);
}























































































bool zut(wfc_blocks_ptr blocks, wfc_blocks_ptr init, uint64_t seed) { return false;}
