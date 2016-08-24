/**
 * Main GCCG program
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "papi.h"
#include "mpi.h"

#include "initialization.h"
#include "compute_solution.h"
#include "finalization.h"

int main(int argc, char *argv[]) {
    int i;

    const int max_iters = 10000;    /// maximum number of iteration to perform

    /** Simulation parameters parsed from the input datasets */
    int nintci, nintcf;    /// internal cells start and end index
    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci, nextcf;
    int **lcc;    /// link cell-to-cell array - stores neighboring information
    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bh, *bl;
    double *bp;    /// Pole coefficient
    double *su;    /// Source values

    double residual_ratio;    /// the ratio between the reference and the current residual
    double *var;    /// the variation vector -> keeps the result in the end

    /** Additional vectors required for the computation */
    double *cgup, *oc, *cnorm;

    /** Command Line arguments **/
    /**gccg <format> <input file> <output prefix>**/
    char *form = argv[1];		// Input file format
    char *file_in = argv[2];	// Input file name
    char *pre = argv[3];		// Output Prefix
 	
    /** PAPI Initialization and Event Creation **/
    long long int values[4];
    int EventSet = PAPI_NULL;

    int Events[4] = {PAPI_L2_TCM, PAPI_L2_TCA, PAPI_L3_TCM, PAPI_L3_TCA};

	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
		printf("Error: PAPI_library_init!");
		exit(1);
	}
	if (PAPI_create_eventset(&EventSet) != PAPI_OK) {
		printf("Error: PAPI_create_eventset");
		exit(1);
	}
	if (PAPI_add_events(EventSet, Events, 4) != PAPI_OK) {
		printf("Error: PAPI_add_events");
		exit(1);
	}

    /********** START INITIALIZATION **********/
    // read-in the input file
    int init_status = initialization(file_in, form, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                     &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &var, &cgup, &oc, 
                                     &cnorm);

    if ( init_status != 0 ) {
        fprintf(stderr, "Failed to initialize data!\n");
        abort();
    } 

   
    /********** END INITIALIZATION **********/

    /** PAPI Start and Verification **/
    if (PAPI_start(EventSet) != PAPI_OK) {
		printf("Error: PAPI_start");
		exit(1);
	}
   
    /********** START COMPUTATIONAL LOOP **********/
    printf("\nComputing Solution...\n");
    int total_iters = compute_solution(max_iters, nintci, nintcf, nextcf, lcc, bp, bs, bw, bl, bn,
                                       be, bh, cnorm, var, su, cgup, &residual_ratio);
    /********** END COMPUTATIONAL LOOP **********/

    /** PAPI End and Verification **/
 	if (PAPI_stop(EventSet, values) != PAPI_OK) {
		printf("Error: PAPI_stop");
		exit(1);
	}
	//ttime = t2 - t1;
    /********** START FINALIZATION **********/
    finalization(file_in, pre, total_iters, residual_ratio, nintci, nintcf, var, cgup, su, lcc, values);
    /********** END FINALIZATION **********/

    free(cnorm);
    free(var);
    free(cgup);
    free(su);
    free(bp);
    free(bh);
    free(bl);
    free(bw);
    free(bn);
    free(be);
    free(bs);

for ( i = 0; i < 6; ++i)
{
   printf("lcc[1][%d] = %d\n", i, lcc[1][i] );
}

   for ( i = nintci; i <= nintcf; i++ ) {
        free(lcc[i]);
    }
   free(lcc);
 
    return 0;
}

