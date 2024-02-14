#include <stdio.h>
#define _gnu_source

#include "bitfield.h"
#include "wfc.h"

#include <omp.h>


bool
solve_cpu(wfc_blocks_ptr blocks)
{
    uint64_t iteration  = 0;
    const uint64_t seed = blocks->seed;
    struct {
        uint32_t gy, x, y, _1;
        uint64_t state;
    } row_changes[blocks->grid_side];

    grd_print(NULL, blocks);
    getchar();

    bool success = false;

    jusqua_la_retraite {
        bool changed = false;

        printf("loc be like   : [gy, gx] [y, x]\n" );
        // choose min entropy 
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
        printf(" choose loc :   [%d, %d] : [%d, %d] = %d\n", gy, gx, y, x, min_entropy );
      
        if( min_entropy != UINT8_MAX ){
            // collapse state
            uint64_t * state = blk_at(blocks, gx, gy, x, y);
            uint64_t collapsed_state = entropy_collapse_state(
                *state, gx, gy, x, y, blocks->seed, iteration);
            *state = collapsed_state;

            changed = grd_propagate_all(blocks, gx, gy, x, y, collapsed_state);
        }
        else{
            success = true;
        }
 

        // 1. collapse
        // 2. propagate
        // 3. check error

        grd_print(NULL, blocks);


        // changed = true;

        iteration += 1;
        getchar();
        if (!changed)
            break;
    }
    // entropy_location loc = min_entropy( blocks );
    entropy_location min_loc;
    min_loc.entropy = UINT8_MAX;
    min_loc.location.x = 0 ; 
    min_loc.location.y = 0;

    // for( uint32_t gy = 0; gy < blocks->grid_side; gy++ ){
    //     for( uint32_t gx = 0; gx < blocks->grid_side; gx++ ){
    //         entropy_location loc = blk_min_entropy(blocks, gx, gy);
    //         if( loc.entropy < min_loc.entropy )
    //             min_loc = loc;
    //     }
    // }
    // uint64_t * state = choose_and_collapse(blocks, 0);

    // *blk_at(ret, gx, gy, x, y) = collapsed;
   

    // entropy_collapse_state(state, )
     

    // for (uint32_t i = 0; i < 3; i++) {
    //     for (uint32_t j = 0; j < 3; j++) {
    //         entropy_location loc = blk_min_entropy(blocks, i, j);
    //         printf("min entropy blk (%u; %u) : [%u; %u] = %u\n",
    //                i, j, loc.location.x, loc.location.y, loc.entropy );
    //     }
    // } 
    grd_print(NULL, blocks);

    return success;
}