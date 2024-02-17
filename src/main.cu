// #define _GNU_SOURCE

#include "bitfield.cuh"
#include "wfc.cuh"
#include "utils.cuh"

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
    
    int max_threads = args.parallel > omp_get_max_threads() ? omp_get_max_threads() : args.parallel;
    
#pragma omp parallel num_threads(max_threads)
    {
        while (!*quit_ptr) {

            wfc_blocks_ptr blocks    = NULL;

            pthread_mutex_lock(&seed_mtx);
            uint64_t next_seed       = 0;
            const bool has_next_seed = try_next_seed(&args.seeds, &next_seed);
            pthread_mutex_unlock(&seed_mtx);

            if (!has_next_seed) {
                __atomic_fetch_or((int*)quit_ptr, true, __ATOMIC_SEQ_CST);
                fprintf(stderr, "T%d : no more seed to try\n", omp_get_thread_num());
                break;
            }

            // wfc_clone_into(&blocks, next_seed, init);
            blocks = cudaCloneToDevice( init );
            
            const bool solved = args.solver(blocks);
            __atomic_add_fetch(iterations_ptr, 1, __ATOMIC_SEQ_CST);

            if (solved && args.output_folder != NULL) {
                __atomic_fetch_or((int*)quit_ptr, true, __ATOMIC_SEQ_CST);
                fputc('\n', stdout);
                wfc_save_into(blocks, args.data_file, args.output_folder);
            }

            else if (solved) {
                __atomic_fetch_or((int*)quit_ptr, 1, __ATOMIC_SEQ_CST);
                fprintf(stdout,"\nT%d : success with result:\n ", omp_get_thread_num());
                // grd_print(NULL, blocks);
            }

            else if (!*quit_ptr) {
                fprintf(stdout, "\rT%d : %.2f%% -> %.2fs\n", omp_get_thread_num(),
                        ((double)(*iterations_ptr) / (double)(max_iterations)) * 100.0,
                        omp_get_wtime() - start);
            }
            // super_safe_free(blocks);
            super_safe_Cudafree(blocks);
            blocks = NULL;
        }
        #pragma omp barrier
    }

    super_safe_free(init);

    return 0;
}
