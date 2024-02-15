#pragma once

#include "types.h"
#include "bitfield.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <strings.h>


static inline wfc_blocks *
super_safe_malloc( uint32_t gs, uint32_t bs)
{

    const uint64_t state_count = gs * gs * bs * bs;
    wfc_blocks *ret = (wfc_blocks *)malloc( sizeof(wfc_blocks) );
    ret->states    = (uint64_t*) malloc( state_count * sizeof(uint64_t) );
    ret->row_masks = (uint64_t*) malloc( gs * bs * sizeof(uint64_t) );
    ret->col_masks = (uint64_t*) malloc( gs * bs * sizeof(uint64_t) );
    ret->blk_masks = (uint64_t*) malloc( gs * gs * sizeof(uint64_t) );
    ret->stack_cells = (vec4*) malloc( (state_count-1) * sizeof(vec4) );

    ret->stack_size = 0;

    uint64_t mask = 0;
    for (uint8_t i = 0; i < bs * bs; i += 1) {
        mask = bitfield_set(mask, i);
    }
    
    for (uint64_t i = 0; i < state_count; i += 1) 
        ret->states[i] = mask;

    for (uint64_t i = 0; i < gs * bs; i += 1) {
        ret->row_masks[i] = mask;
        ret->col_masks[i] = mask;
    }

    for (uint64_t i = 0; i < gs * gs; i += 1) 
        ret->blk_masks[i] = mask;

    return ret; 
}

static inline void
super_safe_free( wfc_blocks * blocks )
{
    free(blocks->stack_cells);
    free(blocks->row_masks);
    free(blocks->col_masks);
    free(blocks->blk_masks);
    free(blocks->states);
    free(blocks);

}
