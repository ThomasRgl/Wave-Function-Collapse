#include <cstdio>
#include <stdint.h>
// #define _GNU_SOURCE

#include "wfc.cuh"
#include "bitfield.cuh"
#include "utils.cuh" 
#include "md5.cuh"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <math.h>

__device__   
uint64_t
entropy_collapse_state(uint64_t state,
                       uint32_t gx, uint32_t gy, uint32_t x, uint32_t y,
                       uint64_t seed,
                       uint64_t iteration)
{
    uint8_t digest[16]     = { 0 };
    uint64_t random_number = 0;
    struct {
        uint32_t gx, gy, x, y;
        uint64_t seed, iteration;
    } random_state = {
        .gx        = gx,
        .gy        = gy,
        .x         = x,
        .y         = y,
        .seed      = seed,
        .iteration = iteration,
    };

    md5((uint8_t *)&random_state, sizeof(random_state), digest);

    uint64_t state_count = bitfield_count(state);
    uint8_t random = (uint8_t)((uint64_t*) &digest)[0];
    uint8_t nth = random%state_count + 1;
    state = bitfield_only_nth_set(state, (uint8_t)nth - 1);

    // printf("new state : %lu (", (uint64_t)log2((double)state)+1);
    // printBinary2(state);
    // printf(")\n");

    return state;
}

__host__ __device__   
uint8_t
entropy_compute(uint64_t state)
{
    return bitfield_count(state);
}

void
wfc_clone_into(wfc_blocks_ptr *const restrict ret_ptr, uint64_t seed, const wfc_blocks_ptr blocks)
{
    const uint32_t gs  = blocks->grid_side;
    const uint32_t bs = blocks->block_side;
    wfc_blocks_ptr ret        = *ret_ptr;

    if (NULL == ret) {
        if (NULL == (ret = super_safe_malloc(gs, bs)) ) {
            fprintf(stderr, "failed to clone blocks structure\n");
            exit(EXIT_FAILURE);
        }
    } else if ( gs != ret->grid_side || bs != ret->block_side) {
        fprintf(stderr, "size mismatch!\n");
        exit(EXIT_FAILURE);
    }
    
    memcpy(ret->states, blocks->states, gs * gs * bs * bs * sizeof(uint64_t) );
    memcpy(ret->row_masks, blocks->row_masks, gs * bs * sizeof(uint64_t) );
    memcpy(ret->col_masks, blocks->col_masks, gs * bs * sizeof(uint64_t) );
    memcpy(ret->blk_masks, blocks->blk_masks, gs * gs * sizeof(uint64_t) );
    ret->seed = seed;
    ret->block_side = (uint8_t)bs;
    ret->grid_side = (uint8_t)gs;

    *ret_ptr       = ret;
}

__device__  
entropy_location
blk_min_entropy(const wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy)
{
    uint8_t min_entropy = UINT8_MAX;
    vec2 the_location   = { 0 };
    uint64_t * blk_state = grd_at(blocks, gx, gy);

    uint32_t bs = blocks->block_side; 
    for (uint32_t y = 0; y < bs; y++) {
        for (uint32_t x = 0; x < bs; x++) {
            uint64_t idx = idx_at(blocks, 0, 0, x, y);
            uint64_t state = blk_state[ idx ] ;
            uint8_t  entropy = entropy_compute(state);
            // printf("idx %d %d = %d == %d\n",y , x, idx, entropy);

            // this condition accept 0
            // if entropy == 0, error will be raise later 
            if( entropy != 1  ){
                if(min_entropy > entropy){
                    min_entropy = entropy;
                    the_location.y = y;
                    the_location.x = x;
                }
            }
        }
    }
    
    entropy_location loc;
    loc.location = the_location;
    loc.entropy = min_entropy;

    // printf("min entropy blk (%u; %u) : [%u; %u] = %u\n",
    //                gy, gx, loc.location.y, loc.location.x, loc.entropy );

    return loc;
}


__device__  
vec4
grd_min_entropy(const wfc_blocks_ptr blocks )
{
    // extern __shared__ entropy_location shared[];
    __shared__ entropy_location shared[8*8];

    shared[threadIdx.y * blockDim.x + threadIdx.x]= blk_min_entropy(blocks, threadIdx.x, threadIdx.y);

    __syncthreads();

    __shared__ vec4 loc;

    if( threadIdx.x == 0 && threadIdx.y == 0 ){
        entropy_location min;
        min.entropy = UINT8_MAX;
        min.location.x = UINT8_MAX;
        min.location.y = UINT8_MAX;
        uint32_t gx;
        uint32_t gy;
        for (uint32_t i = 0; i < blockDim.x * blockDim.y; i++) {
            if( shared[i].entropy < min.entropy ){
                min = shared[i];
                gx = i % blockDim.x;
                gy = i / blockDim.x;
            }
        }
        loc.gx = gx;
        loc.gy = gy;
        loc.y = min.location.y;
        loc.x = min.location.x;
    }

    __syncthreads();

    // printf("T%u : min entropy grd : [%u; %u] [%u; %u] = %u\n", threadIdx.x + threadIdx.y * blockDim.x,
    // loc.gy, loc.gx, loc.y, loc.x, min.entropy );

    return loc;
}


__device__  
bool
blk_propagate(wfc_blocks_ptr blocks,
              uint32_t gx, uint32_t gy,
              uint64_t collapsed)
{
    uint32_t x = threadIdx.x;
    uint32_t y = threadIdx.y;
    
    
    if ( (blocks->blk_masks[gy * blocks->block_side + gx] & collapsed) == 0) {
        // if (threadIdx.x == 0 && threadIdx.y == 0) {
        //     printf("error in (mask) block propagation in block (%u, %u)\n", gy, gx);
        // }
        return true;
    }

    blocks->blk_masks[gy * blocks->block_side + gx] &= ~collapsed;

    uint64_t idx = idx_at(blocks, gx, gy, x, y);
    uint8_t entropy = entropy_compute(blocks->states[idx]);
    // if the cell is not collapsed yet
    if (entropy > 1) {
        blocks->states[idx] &= ~(collapsed);
        // if the new entropy is 1 add the cell to the stack to propagate it later
        entropy = entropy_compute(blocks->states[idx]);
        if (entropy == 1) {
            vec4 coord = {gx, gy, x, y};
            uint32_t unique_stack_id = atomicAdd(&blocks->stack_size, 1);
            blocks->stack_cells[unique_stack_id] = coord;
            // printf("T%u : add stack (%d)  (blk)\n", threadIdx.x + threadIdx.y * blockDim.x, blocks->stack_size);
        }
        else if (entropy == 0) {
            // fprintf(stderr, "error in (mask) block propagation in block (%u, %u)\n", gy, gx);
            return true;
        }
    }

    return false;
}

__device__  
bool
grd_propagate_column(wfc_blocks_ptr blocks,
                  uint32_t gx, uint32_t __, uint32_t x, uint32_t _,
                  uint64_t collapsed)
{
    uint32_t gy = threadIdx.x;
    uint32_t  y = threadIdx.y;

    if ( (blocks->col_masks[gx * blocks->block_side + x] & collapsed) == 0) {
        // if (threadIdx.x == 0 && threadIdx.y == 0) {
        //     printf("error in (mask) column propagation in column block %u in column %u\n", gx, x);
        // }
        return true;
    }

    blocks->col_masks[gx * blocks->block_side + x] &= ~collapsed;

    uint64_t idx = idx_at(blocks, gx, gy, x, y);
    uint8_t entropy = entropy_compute(blocks->states[idx]);
    // if the cell is not collapsed yet
    if (entropy > 1) {
        blocks->states[idx] &= ~(collapsed);
        // if the new entropy is 1 add the cell to the stack to propagate it later
        entropy = entropy_compute(blocks->states[idx]);
        if (entropy == 1) {
            vec4 coord = {gx, gy, x, y};
            uint32_t unique_stack_id = atomicAdd(&blocks->stack_size, 1);
            blocks->stack_cells[unique_stack_id] = coord;

            // printf("T%u : add stack (%d)  (col)\n", threadIdx.x + threadIdx.y * blockDim.x, blocks->stack_size);
        }
        else if (entropy == 0) {
            // fprintf(stderr, "error in (mask) col propagation in block (%u, %u)\n", gy, gx);
            return true;
        }
    }

    return false;
}


__device__  
bool
grd_propagate_row(wfc_blocks_ptr blocks, uint32_t __, uint32_t gy,
                     uint32_t _, uint32_t y, uint64_t collapsed)
{
    uint32_t gx = threadIdx.x;
    uint32_t  x = threadIdx.y;

    if ( (blocks->row_masks[gy * blocks->block_side + y] & collapsed) == 0) {
        // if (threadIdx.x == 0 && threadIdx.y == 0) {
        //     printf("error in (mask) row propagation in row block %u in row %u\n", gy, y);
        // }
        return true;
    }

    blocks->row_masks[gy * blocks->block_side + y] &= ~collapsed;

    uint64_t idx = idx_at(blocks, gx, gy, x, y);
    uint8_t entropy = entropy_compute(blocks->states[idx]);
    // if the cell is not collapsed yet
    if (entropy > 1) {
        blocks->states[idx] &= ~(collapsed);
        // if the new entropy is 1 add the cell to the stack to propagate it later
        entropy = entropy_compute(blocks->states[idx]);
        if (entropy == 1) {
            vec4 coord = {gx, gy, x, y};
            uint32_t unique_stack_id = atomicAdd(&blocks->stack_size, 1);
            blocks->stack_cells[unique_stack_id] = coord;

            // printf("T%u : add stack (%d)  (row)\n", threadIdx.x + threadIdx.y * blockDim.x, blocks->stack_size);
        }
        else if (entropy == 0) {
            // fprintf(stderr, "error in (mask) row propagation in block (%u, %u)\n", gy, gx);
            return true;
        }
    }

    return false;
}

__device__  
bool
grd_propagate_all(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y, uint64_t collapsed)
{
    __shared__ bool shared_error[8*8];
    bool error = false;

    if( threadIdx.x == 0 && threadIdx.y == 0 )
        blocks->stack_size = 0;
    __syncthreads(); // Attendre de l'ecriture de stack size


    // printf("T%u : collapsed (threads): %lu () at : [%u, %u] [%u, %u]\n", threadIdx.x + threadIdx.y * blockDim.x, (uint64_t)log2((double)collapsed)+1, gy, gx, y, x);
    
    // propagate the initial cell
    error |= grd_propagate_column(blocks, gx, gy, x, y, collapsed);
    error |= grd_propagate_row(blocks, gx, gy, x, y, collapsed);
    __syncthreads(); 
    error |= blk_propagate(blocks, gx, gy, collapsed);

    shared_error[ threadIdx.x + threadIdx.y * blockDim.x ] = error;
    __syncthreads();// 
        
    for(uint32_t i = 0; i < blockDim.x * blockDim.y; i++){
        error |= shared_error[i];        
        // printf("T%u : [%d]  %lu\n", threadIdx.x + threadIdx.y * blockDim.x, i, shared_error[i]);
    }
    
    while ( !error && blocks->stack_size > 0) {

        // if( threadIdx.x == 0 && threadIdx.y == 0 ){
        //     printf("stack (%d): ", blocks->stack_size);
        //     for (uint32_t i = 0; i < blocks->stack_size; i++) {
        //         printf("([%u, %u] [%u, %u]) ", blocks->stack_cells[i].gy, blocks->stack_cells[i].gx, blocks->stack_cells[i].y, blocks->stack_cells[i].x);
        //     }
        //     printf("\n");
        // }

        // pop the cell from the stack
        vec4 coord = blocks->stack_cells[0];
        x = coord.x;
        y = coord.y;
        gx = coord.gx;
        gy = coord.gy;

        __syncthreads(); // On attend que tout les threads lises stack[0]
        if( threadIdx.x == 0 && threadIdx.y == 0 ){
            for (uint32_t i = 0; i < blocks->stack_size; i++) {
                blocks->stack_cells[i] = blocks->stack_cells[i + 1];
            }
            blocks->stack_size--;
        }

        // get the new collapsed state
        collapsed = *blk_at(blocks, gx, gy, x, y);

        // if( threadIdx.x == 0 && threadIdx.y == 0 ){
        //     printf("collapsed (stack): %lu (", (uint64_t)__ffsll((long long)collapsed) );
        //     printBinary2(collapsed);
        //     printf(") at : [%u, %u] [%u, %u]\n", gy, gx, y, x);
        // }

        // printf("T%u : collapsed (threads): at : [%u, %u] [%u, %u]\n", threadIdx.x + threadIdx.y * blockDim.x, gy, gx, y, x);

        // __syncthreads();
        // propagate the new cell
        error |= grd_propagate_column(blocks, gx, gy, x, y, collapsed);
        error |= grd_propagate_row(blocks, gx, gy, x, y, collapsed);
        __syncthreads();
        error |= blk_propagate(blocks, gx, gy, collapsed);

        shared_error[ threadIdx.x + threadIdx.y * blockDim.x ] = error;
        __syncthreads(); // attendre fin shared error ecriture


        if( threadIdx.x == 0 && threadIdx.y == 0 ){
            *blk_at(blocks, gx, gy, x, y) = collapsed; // Useless ????
        }
        __syncthreads(); // wait useless
            
        for(uint32_t i = 0; i < blockDim.x * blockDim.y; i++){
            error |= shared_error[i];        
        }
        // __syncthreads(); // wait 
    }

    return error;
}

__host__ __device__  
void printBinary2(uint64_t number) {
    // Determine the number of bits in u64_t
    int numBits = sizeof(uint64_t) + 1;

    // Loop through each bit in the number, starting from the most significant bit
    for (int i = numBits - 1; i >= 0; i--) {
        // Use a bitwise AND operation to check the value of the current bit
        if ((number & (1ULL << i)) != 0) {
            printf("1");
        } else {
            printf("0");
        }

        // Add space for better readability
        // if (i % 8 == 0) {
        //     printf(" ");
        // }
    }
    // printf("\n");
}

__host__ __device__  
void grd_print(FILE *const file, const wfc_blocks_ptr block){
    
    uint8_t gs = block->grid_side;
    uint8_t bs = block->block_side;

    for(uint32_t ii = 0; ii < gs; ii++){

        for(uint32_t i = 0; i < gs; i++){
            printf( "+");
            for(uint32_t j = 0; j < bs; j++){
                printf( "----------+");
            }
            printf( "   ");
        }
        printf( "\n");

        for(uint32_t jj = 0; jj < bs; jj++){

            for(uint32_t i = 0; i < gs; i++){
                printf( "|");
                for(uint32_t j = 0; j < bs; j++){
                    const uint64_t collapsed = *blk_at(block, i,ii,j,jj);
                    if( bitfield_count(collapsed) == 1){

                        printf( "    %2lu   ", 
                        #ifdef __CUDA_ARCH__
                        (uint64_t)__ffsll((long long)collapsed));
                        #else
                        (uint64_t)log2((double)collapsed)+1 );
                        #endif                    
                        // printBinary2(collapsed);
                    }
                    else
                        printBinary2(collapsed);
                    // printf( "   %lu     ", idx_at(block, i,ii,j,jj) );
                    printf( " |");
                    // printf( "%3lu|");
                }
                printf( "   ");
            }
            printf( "\n");

            for(int i = 0; i < gs; i++){
                printf( "+");
                for(int j = 0; j < bs; j++){
                    printf( "----------+");
                }
                printf( "   ");
            }
            printf( "\n");
        }
        printf( "\n");
    }

}





/////////////////// HOST PROPAGE ////////////////////

__host__  
bool
host_blk_propagate(wfc_blocks_ptr blocks,
              uint32_t gx, uint32_t gy,
              uint64_t collapsed)
{
    
    
    if ( (blocks->blk_masks[gy * blocks->block_side + gx] & collapsed) == 0) {
        // fprintf(stderr, "error in (mask) block propagation in block (%u, %u)\n", gy, gx);
        return true;
    }

    blocks->blk_masks[gy * blocks->block_side + gx] &= ~collapsed;

    for (uint32_t y = 0; y < blocks->block_side; y++) {
        for (uint32_t x = 0; x < blocks->block_side; x++) {

            uint64_t idx = idx_at(blocks, gx, gy, x, y);
            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    vec4 coord = {gx, gy, x, y};
                    blocks->stack_cells[blocks->stack_size] = coord;
                    blocks->stack_size++;
                }
                else if (entropy == 0) {
                    // fprintf(stderr, "error in (mask) block propagation in block (%u, %u)\n", gy, gx);
                    return true;
                }
            }
        }
    }

    return false;
}

__host__  
bool
host_grd_propagate_column(wfc_blocks_ptr blocks,
                  uint32_t gx, uint32_t __, uint32_t x, uint32_t _,
                  uint64_t collapsed)
{

    if ( (blocks->col_masks[gx * blocks->block_side + x] & collapsed) == 0) {
        // fprintf(stderr, "error in (mask) column propagation in column block %u in column %u\n", gx, x);
        return true;
    }

    blocks->col_masks[gx * blocks->block_side + x] &= ~collapsed;

    for (uint32_t gy = 0; gy < blocks->grid_side; gy++) {
        for (uint32_t y = 0; y < blocks->block_side; y++) {
            uint64_t idx = idx_at(blocks, gx, gy, x, y);

            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    vec4 coord = {gx, gy, x, y};
                    blocks->stack_cells[blocks->stack_size] = coord;
                    blocks->stack_size++;
                }
                else if (entropy == 0) {
                    // fprintf(stderr, "error in (mask) col propagation in block (%u, %u)\n", gy, gx);
                    return true;
                }
                
            }
        }
    }

    return false;
}


__host__  
bool
host_grd_propagate_row(wfc_blocks_ptr blocks, uint32_t __, uint32_t gy,
                     uint32_t _, uint32_t y, uint64_t collapsed)
{

    if ( (blocks->row_masks[gy * blocks->block_side + y] & collapsed) == 0) {
        // fprintf(stderr, "error in (mask) row propagation in row block %u in row %u\n", gy, y);
        return true;
    }

    blocks->row_masks[gy * blocks->block_side + y] &= ~collapsed;

    for (uint32_t gx = 0; gx < blocks->grid_side; gx++) {
        for (uint32_t x = 0; x < blocks->block_side; x++) {
            // uint64_t idx = i * blocks->grid_side * blk_size + gy * blk_size + j * blocks->block_side + y;
            uint64_t idx = idx_at(blocks, gx, gy, x, y);
            uint8_t entropy = entropy_compute(blocks->states[idx]);
            // if the cell is not collapsed yet
            if (entropy > 1) {
                blocks->states[idx] &= ~(collapsed);
                // if the new entropy is 1 add the cell to the stack to propagate it later
                entropy = entropy_compute(blocks->states[idx]);
                if (entropy == 1) {
                    vec4 coord = {gx, gy, x, y};
                    blocks->stack_cells[blocks->stack_size] = coord;
                    blocks->stack_size++;
                }
                else if (entropy == 0) {
                    // fprintf(stderr, "error in (mask) row propagation in block (%u, %u)\n", gy, gx);
                    return true;
                }
            }
        }
    }

    return false;
}

__host__  
bool
host_grd_propagate_all(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y, uint64_t collapsed)
{
    bool error = false;

    // propagate the initial cell
    error |= host_blk_propagate(blocks, gx, gy, collapsed);
    error |= host_grd_propagate_column(blocks, gx, gy, x, y, collapsed);
    error |= host_grd_propagate_row(blocks, gx, gy, x, y, collapsed);
    *blk_at(blocks, gx, gy, x, y) = collapsed;
    
    while ( !error && blocks->stack_size > 0) {

        // printf("stack (%d): ", blocks->stack_size);
        // for (int i = 0; i < blocks->stack_size; i++) {
        //     printf("%d ", blocks->stack_cells[i]);
        // }
        // printf("\n");

        // pop the cell from the stack
        vec4 coord = blocks->stack_cells[0];
        x = coord.x;
        y = coord.y;
        gx = coord.gx;
        gy = coord.gy;
        for (uint32_t i = 0; i < blocks->stack_size; i++) {
            blocks->stack_cells[i] = blocks->stack_cells[i + 1];
        }
        blocks->stack_size--;

        // get the new collapsed state
        collapsed = *blk_at(blocks, gx, gy, x, y);
        // if(collapsed == 0)
        //     return error;

        // printf("collapsed (stack): %lu (", collapsed);
        // printBinary2(collapsed);
        // printf(") at : %lu, %lu, %lu, %lu\n", gy, gx, y, x);

        // propagate the new cell
        error |= host_blk_propagate(blocks, gx, gy, collapsed);
        error |= host_grd_propagate_column(blocks, gx, gy, x, y, collapsed);
        error |= host_grd_propagate_row(blocks, gx, gy, x, y, collapsed);
        *blk_at(blocks, gx, gy, x, y) = collapsed;
    }

    return error;
}

__host__  
bool
verify_block(wfc_blocks_ptr blocks )
{
    uint32_t gs = blocks->grid_side;
    uint32_t bs = blocks->block_side;

    for (uint32_t gy = 0; gy < gs; gy++) {
        for (uint32_t gx = 0; gx < gs; gx++) {

            uint64_t mask = 0;
            for (uint32_t y = 0; y < bs; y++) {
                for (uint32_t x = 0; x < bs; x++) {
                    uint32_t idx = (uint32_t)idx_at(blocks, gx, gy, x, y);
                    uint64_t state = blocks->states[idx];
                    if( (mask & state) != 0){
                        printf("error block at [%u %u] [%u %u] = %lu\n", gy, gx, y,x,  (uint64_t)log2((double)state)+1); 
                        exit(-1);
                    }
                    mask |= blocks->states[idx];
                }
            }
            
        }
    }

    for (uint32_t gy = 0; gy < gs; gy++) {
        for (uint32_t y = 0; y < bs; y++) {

            uint64_t mask = 0;
            for (uint32_t gx = 0; gx < gs; gx++) {
                for (uint32_t x = 0; x < bs; x++) {
                    uint32_t idx = (uint32_t)idx_at(blocks, gx, gy, x, y);
                    uint64_t state = blocks->states[idx];
                    if( (mask & state) != 0){
                        printf("error COL at [%u %u] [%u %u] = %lu\n", gy, gx, y,x, (uint64_t)log2((double)state)+1 ); 
                        exit(-1);
                    }
                    mask |= blocks->states[idx];
                }
            }
            
        }
    }

    for (uint32_t gx = 0; gx < gs; gx++) {
        for (uint32_t x = 0; x < bs; x++) {

            uint64_t mask = 0;
            for (uint32_t gy = 0; gy < gs; gy++) {
                for (uint32_t y = 0; y < bs; y++) {
                    uint32_t idx = (uint32_t)idx_at(blocks, gx, gy, x, y);
                    uint64_t state = blocks->states[idx];
                    if( (mask & state) != 0){
                        printf("error COL at [%u %u] [%u %u] = %lu\n", gy, gx, y,x, (uint64_t)log2((double)state)+1 ); 
                        exit(-1);
                    }
                    mask |= blocks->states[idx];
                }
            }
            
        }
    }
    

    return true;
}
