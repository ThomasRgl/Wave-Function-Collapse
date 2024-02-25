// #define _GNU_SOURCE

#include "bitfield.cuh"
#include "wfc.cuh"
#include "utils.cuh"

#include <cstdio>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <pthread.h>

int
main(int argc, char **argv)
{
    omp_set_dynamic(false);

    wfc_args args             = wfc_parse_args(argc, argv);
    const wfc_blocks_ptr init = wfc_load(0, args.data_file);

    bool quit                = false;
    uint64_t iterations      = 0;
    pthread_mutex_t seed_mtx = PTHREAD_MUTEX_INITIALIZER;

    bool *volatile const quit_ptr           = &quit;
    uint64_t *volatile const iterations_ptr = &iterations;

    const uint64_t max_iterations = count_seeds(args.seeds);
    const double start            = omp_get_wtime();

    int nb_device;
    cudaGetDeviceCount(&nb_device);
    printf("%d devices founds\n", nb_device);
    
    // int max_threads = 1;
    uint64_t max_threads = 4 > nb_device ? (uint64_t)nb_device : (uint64_t)4;
    if( omp_get_max_threads() < (int)max_threads ){
        printf(" not enough OMP threads available\n");
        return 0;
    }

    printf("running on %lu devices \n", max_threads);

    wfc_blocks_ptr * d_blocks_list = (wfc_blocks_ptr *) malloc( max_threads * sizeof(wfc_blocks_ptr));
    wfc_blocks_ptr * d_init_list   = (wfc_blocks_ptr *) malloc( max_threads * sizeof(wfc_blocks_ptr));

    // const wfc_blocks_ptr d_init = wfc_clone_HTD( init);

    for (int i = 0; i < max_threads; i++) { 
        cudaSetDevice((int)i);
        // printf("SET DEVIDE %d\n", (int)i);
        // blocks_list[i] = cloneToDevice( init, 0 );
        d_init_list[i]   = wfc_clone_HTD(init);
        d_blocks_list[i] = wfc_clone_HTD(init);
    }
    
    
#pragma omp parallel num_threads(max_threads)
    {
        while (!*quit_ptr) {

            cudaSetDevice(omp_get_thread_num());
            // printf("DEVICE : %d\n", omp_get_thread_num() );
            wfc_blocks_ptr d_blocks = d_blocks_list[omp_get_thread_num()];
            wfc_blocks_ptr d_init   = d_init_list[omp_get_thread_num()];

            // pthread_mutex_lock(&seed_mtx);
            // uint64_t next_seed       = 0;
            // const bool has_next_seed = try_next_seed(&args.seeds, &next_seed);
            // pthread_mutex_unlock(&seed_mtx);

            pthread_mutex_lock(&seed_mtx);
            uint64_t seed_list[args.parallel];
            bool has_next_seed = true;
            int nb_seed = 0;

            // tant qu'il reste des seeds a explorer et que i << args.par
            for( int i = 0; i < args.parallel && has_next_seed; i++) {
                has_next_seed = try_next_seed(&args.seeds, &seed_list[i]);
                nb_seed ++;
            }

            // si on a depassé le nombre de seed existante
            if(!has_next_seed){
                nb_seed--;
            }

            pthread_mutex_unlock(&seed_mtx);

            // si il ya des seed a explorer
            if (nb_seed == 0) {
                __atomic_fetch_or((int*)quit_ptr, true, __ATOMIC_SEQ_CST);
                fprintf(stderr, "T%d : no more seed to try\n", omp_get_thread_num());
                break;
            }

            // wfc_clone_into(&blocks, next_seed, init);
            // blocks = cudaCloneToDevice( init, next_seed );
            
            wfc_blocks_ptr h_blocks = solve_cuda(d_blocks, d_init, 
                    seed_list, init->grid_side, init->block_side, nb_seed);

            bool solved = false;

            if (h_blocks != NULL) {
                solved = h_blocks->solved;
            }

            // __atomic_add_fetch(iterations_ptr, 1, __ATOMIC_SEQ_CST);
            __atomic_add_fetch(iterations_ptr, nb_seed, __ATOMIC_SEQ_CST);


            if (solved && args.output_folder != NULL) {
                __atomic_fetch_or((int*)quit_ptr, true, __ATOMIC_SEQ_CST);
                fputc('\n', stdout);
                #pragma omp single nowait
                {
                    wfc_save_into(h_blocks, args.data_file, args.output_folder);
                }
            }

            else if (solved) {
                __atomic_fetch_or((int*)quit_ptr, 1, __ATOMIC_SEQ_CST);
                fprintf(stdout,"\nT%d : success with result:\n ", omp_get_thread_num());
                // grd_print(NULL, blocks);
            }

            if (h_blocks != NULL) {
                super_safe_free(h_blocks);
            }

            //else if (!*quit_ptr) {
                fprintf(stdout, "\rT%d : %.2f%% -> %.2fs\n", omp_get_thread_num(),
                        ((double)(*iterations_ptr) / (double)(max_iterations)) * 100.0,
                        omp_get_wtime() - start);
            //}
            // super_safe_free(blocks);
            // d_blocks = NULL;
        }
        #pragma omp barrier
    }

    for (uint64_t i = 0; i < max_threads; i++){
        cudaSetDevice((int)i);
        super_safe_Cudafree(d_blocks_list[i]);
        super_safe_Cudafree(d_init_list[i]);
    }
    free(d_blocks_list);
    free(d_init_list);

    super_safe_free(init);

    return 0;
}