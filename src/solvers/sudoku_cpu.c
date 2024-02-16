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
     
    grd_print(NULL, blocks);
    getchar();

    return success;
}
