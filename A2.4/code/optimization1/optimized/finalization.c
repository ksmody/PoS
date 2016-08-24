/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "util_write_files.h"

void finalization(char* file_in, int nprocs, int myrank, int total_iters,
		double residual_ratio, int nintci, int nintcf, double* var,
		int* local_global_index, int* global_local_index) {

	int i, nintcf_global, *lg_p, *rcounts, *rdispls;
	double *var_1, *var_2;

	/* Summing up the total Internal Cells in Domain in Rank 0*/
	nintcf++;
	MPI_Reduce(&nintcf, &nintcf_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Setting up the variables for the combining of VAR array*/
	if (myrank == 0) {
		var_1 = malloc(sizeof(double) * nintcf_global);
		var_2 = malloc(sizeof(double) * nintcf_global);
		lg_p = malloc(sizeof(int) * nintcf_global);
		rcounts = malloc(sizeof(double) * nprocs);
		rdispls = malloc(sizeof(double) * nprocs);
	}

	MPI_Gather(&nintcf, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (myrank == 0) {
		rdispls[0] = 0;
		for (i = 1; i < nprocs; ++i) {
			rdispls[i] = rdispls[i - 1] + rcounts[i - 1];
		}

	}

	/*Gathering all necessary data for forming Global VAR*/
	MPI_Gatherv(local_global_index, nintcf, MPI_INT, lg_p, rcounts, rdispls,
			MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(var, nintcf, MPI_DOUBLE, var_1, rcounts, rdispls, MPI_DOUBLE, 0,
			MPI_COMM_WORLD);

	if (myrank == 0) {
		/*Re-Ordering the VAR values based on local to global index*/
		for (i = 0; i < nintcf_global; ++i) {
			var_2[lg_p[i]] = var_1[i];
		}

		nintcf_global--;

		/*Creating the Summary Output File*/
		char file_out[100];
		sprintf(file_out, "%s_summary.out", file_in);

		/*Function call for displaying Output and .out file*/
		int status = store_simulation_stats(file_in, file_out, nintci,
				nintcf_global, var_2, total_iters, residual_ratio);

		if (status != 0)
			fprintf(stderr, "Error when trying to write to file %s\n",
					file_out);

	}
}
