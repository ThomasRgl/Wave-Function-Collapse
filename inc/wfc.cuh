#pragma once

#include "types.cuh"

#include <stdbool.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

#define jusqua_la_retraite for (;;)

/// Parses the arguments, prints the help message if needed and abort on error.
wfc_args wfc_parse_args(int argc, char **argv);

/// Get the next seed to try. If there are no more seeds to try, it will exit the process.
bool try_next_seed(seeds_list *restrict *const, uint64_t *restrict);

/// Count the total number of seeds.
uint64_t count_seeds(const seeds_list *restrict const);

/// Load the positions from a file. You must free the thing yourself. On error
/// kill the program.
wfc_blocks_ptr wfc_load(uint64_t, const char *);

/// Clone the blocks structure. You must free the return yourself.
void wfc_clone_into(wfc_blocks_ptr *const restrict, uint64_t, const wfc_blocks_ptr);

/// Save the grid to a folder by creating a new file or overwrite it, on error kills the program.
void wfc_save_into(const wfc_blocks_ptr, const char data[], const char folder[]);

__host__ __device__  static inline uint64_t
wfc_control_states_count(uint64_t grid_size, uint64_t block_size)
{
    return grid_size * grid_size * block_size * block_size;
}

__host__ __device__ static inline uint64_t 
idx_at(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
{
    uint64_t idx = gy * blocks->grid_side * blocks->block_side * blocks->block_side +
                   gx * blocks->block_side * blocks->block_side + y * blocks->block_side + x;
    return idx;
}

__host__ __device__ static inline uint64_t *
grd_at(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy)
{
    uint64_t idx = idx_at(blocks, gx, gy, 0, 0);
    return &blocks->states[idx];
}

__host__ __device__ static inline uint64_t *
blk_at(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y)
{
    uint64_t idx = idx_at(blocks, gx, gy, x, y);

    return &blocks->states[idx];
}





__host__ __device__ void printBinary2(uint64_t number) ;

// void blk_print(FILE *const, const wfc_blocks_ptr block, uint32_t gx, uint32_t gy);
__host__ __device__ void grd_print(FILE *const, const wfc_blocks_ptr block);

// Entropy functionsvec4
__device__ vec4 grd_min_entropy(const wfc_blocks_ptr blocks );
__device__ entropy_location blk_min_entropy(const wfc_blocks_ptr block, uint32_t gx, uint32_t gy);
__host__ __device__ uint8_t entropy_compute(uint64_t);
__device__ uint64_t entropy_collapse_state(uint64_t, uint32_t, uint32_t, uint32_t, uint32_t, uint64_t, uint64_t);

// Propagation functions

__host__ bool host_grd_propagate_all(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y, uint64_t collapsed);
__host__ bool host_blk_propagate(wfc_blocks_ptr, uint32_t, uint32_t, uint64_t);
__host__ bool host_grd_propagate_column(wfc_blocks_ptr, uint32_t, uint32_t, uint32_t, uint32_t, uint64_t);
__host__ bool host_grd_propagate_row(wfc_blocks_ptr, uint32_t, uint32_t, uint32_t, uint32_t, uint64_t);
__host__ bool verify_block(wfc_blocks_ptr blocks );

__device__  bool grd_propagate_all(wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy, uint32_t x, uint32_t y, uint64_t collapsed);
__device__  bool blk_propagate(wfc_blocks_ptr, uint32_t, uint32_t, uint64_t);
__device__  bool grd_propagate_column(wfc_blocks_ptr, uint32_t, uint32_t, uint32_t, uint32_t, uint64_t);
__device__  bool grd_propagate_row(wfc_blocks_ptr, uint32_t, uint32_t, uint32_t, uint32_t, uint64_t);



// Solvers
// bool solve_cpu(wfc_blocks_ptr);
// bool solve_openmp(wfc_blocks_ptr);
// bool solve_target(wfc_blocks_ptr);
wfc_blocks_ptr solve_cuda(wfc_blocks_ptr d_blocks, wfc_blocks_ptr d_init, uint64_t * seed_list, uint32_t gs, uint32_t bs, uint64_t p);
bool zut(wfc_blocks_ptr blocks, wfc_blocks_ptr init, uint64_t seed);



static const wfc_solver solvers[] = {
    // { "cpu", solve_cpu },
    // { "omp", solve_openmp },
    // { "target", solve_target },
    { "cuda", zut }
};
