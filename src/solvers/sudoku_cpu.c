#include <stdio.h>
#define _gnu_source

#include "bitfield.h"
#include "wfc.h"

#include <omp.h>

entropy_location
blk_min_entropy(const wfc_blocks_ptr blocks, uint32_t gx, uint32_t gy)
{
    uint64_t minima;
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
    
    // uint32_t candidate = bitfield_count(minima);
    entropy_location loc;
    loc.location = the_location;
    loc.entropy = min_entropy;

    // printf("min entropy blk (%u; %u) : [%u; %u] = %u\n",
    //                gy, gx, loc.location.y, loc.location.x, loc.entropy );

    return loc;
}


vec4
grd_min_entropy(const wfc_blocks_ptr blocks )
{
    uint8_t min_entropy = UINT8_MAX;
    uint32_t x = 0 ; 
    uint32_t y = 0;
    uint32_t gy = 0;
    uint32_t gx = 0;
    for( uint32_t gy_ = 0; gy_ < blocks->grid_side; gy_++ ){
        for( uint32_t gx_ = 0; gx_ < blocks->grid_side; gx_++ ){
            entropy_location loc = blk_min_entropy(blocks, gx_, gy_);
            if( loc.entropy < min_entropy ){
                min_entropy = loc.entropy;
                x = loc.location.x;
                y = loc.location.y;
                gy = gy_;
                gx = gx_;
            }
        }
    }

    // printf("choosed min entropy (%u; %u) : [%u; %u] = %u ",
    //                gy, gx, y, x, min_entropy );
    // printBinary2(min_entropy);
    // printf("\n");



    vec4 loc;
    loc.gx = gx;
    loc.gy = gy;
    loc.y = y;
    loc.x = x;

    return loc;
}


bool
solve_cpu(wfc_blocks_ptr blocks)
{
    uint64_t iteration  = 0;
    const uint64_t seed = blocks->seed;
    struct {
        uint32_t gy, x, y, _1;
        uint64_t state;
    } row_changes[blocks->grid_side];

    // grd_print(NULL, blocks);
    // getchar();

    bool success = false;

    jusqua_la_retraite {

        // grd_print(NULL, blocks);

        vec4 loc = grd_min_entropy(blocks);
        uint64_t * state = blk_at(blocks, loc.gx, loc.gy,
                                      loc.x, loc.y);
        // printf(" choose loc :   [%d, %d] : [%d, %d] = %lu\n", 
        // loc.gy, loc.gx, loc.y, loc.x, *state );
      
        if( *state == UINT8_MAX ){
            success = true;
            // printf("success\n");
            break;
        }

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
     
    // grd_print(NULL, blocks);
    // getchar();

    return success;
}
