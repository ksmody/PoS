/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util_read_files.h"
#include "util_write_files.h"
#include "initialization.h"

#include "test_functions.h"
#include "metis.h"
#include "mpi.h"
#include "papi.h"

#define MIN(a,b) (((a)<(b))?(a):(b))

int initialization(char* file_in, char* part_type, char* read_type, int nprocs,
		int myrank, int* nintci, int* nintcf, int* nextci, int* nextcf,
		int*** lcc, double** bs, double** be, double** bn, double** bw,
		double** bl, double** bh, double** bp, double** su, int* points_count,
		int*** points, int** elems, double** var, double** cgup, double** oc,
		double** cnorm, int** local_global_index) {

	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
		printf("Error: PAPI_library_init!");
		exit(1);
	}
	long long int t, t1, t2;
	t1 = PAPI_get_virt_usec();

	int local_intc, local_extc;
	int input_key = 1 , part_key = 1 , read_key = 1;
	/********** START INITIALIZATION **********/
	int *nintcf_a, *displs, *nintcf_a2, *displs2;
	int **lcc_m;
	double *bs_m, *be_m, *bn_m, *bw_m, *bl_m, *bh_m, *bp_m, *su_m;
	int *lcc_a_m, *lcc_a, *nextcf_a;
	int* l_g;
	int i, j;

	if (strcmp(file_in, "tjunc") == 0)
		input_key = 1;
	else if (strcmp(file_in, "drall") == 0)
		input_key = 2;
	else if (strcmp(file_in, "pent") == 0)
		input_key = 3;
	else if (strcmp(file_in, "cojack") == 0)
		input_key = 4;

	if (strcmp(read_type, "oneread") == 0) {
		read_key = 1;

		if (myrank == 0) {

			double *bp_m0, *bs_m0, *be_m0, *bn_m0, *bw_m0, *bl_m0, *bh_m0,
					*su_m0;
			int nintcf_m, nextci_m, nextcf_m, lcc_ij;
			int **lcc_m, **points_l;
			int *epart_ghost;
			int *g_l, *g_l_ghost;

			int r, q, s = 1, k = 0;

			int f_status = read_binary_geo(file_in, nintci, &nintcf_m,
					&nextci_m, &nextcf_m, &lcc_m, &bs_m0, &be_m0, &bn_m0,
					&bw_m0, &bl_m0, &bh_m0, &bp_m0, &su_m0, points_count,
					&points_l, elems);
			if (f_status != 0)
				return f_status;

			g_l = (int*) malloc(sizeof(int) * (nintcf_m + 1));
			g_l_ghost = (int*) malloc(sizeof(int) * (nextcf_m + 1));
			epart_ghost = (int*) malloc(sizeof(int) * (nextcf_m + 1));

			lcc_a_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
			l_g = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);

			bp_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			bs_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			be_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			bn_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			bw_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			bl_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			bh_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
			su_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));

			displs = (int*) malloc(sizeof(int) * nprocs);
			displs2 = (int*) malloc(sizeof(int) * nprocs);
			nintcf_a = (int*) calloc(sizeof(int), nprocs);
			nintcf_a2 = (int*) calloc(sizeof(int), nprocs);
			nextcf_a = (int*) calloc(sizeof(int), nprocs);

			idx_t ne = nintcf_m + 1;
			idx_t nn = *points_count;
			idx_t ncommon = 4;
			idx_t nparts = nprocs;

			idx_t *eptr = (idx_t*) malloc((ne + 1) * sizeof(idx_t*));
			idx_t *eind = (idx_t*) malloc(ne * 8 * sizeof(idx_t*));
			idx_t *objval = (idx_t*) malloc(1 * sizeof(idx_t*));
			idx_t *epart = (idx_t*) malloc((ne) * sizeof(idx_t*));
			idx_t *npart = (idx_t*) malloc((nn) * sizeof(idx_t*));

			for (i = 0; i < (ne + 1); i++) {
				eptr[i] = i * 8;
			}

			for (i = 0; i < (ne * 8); i++) {
				eind[i] = (*elems)[i];
			}

			/**********METIS DISTRIBUTION*************/
			if (strcmp(part_type, "dual") == 0) {
				part_key = 2;
				if (METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL,
						&ncommon, &nparts, NULL, NULL, objval, epart, npart)
						!= METIS_OK)
					printf("ERROR IN METIS DUAL DISTRIBUTION\n");
				else
					printf("METIS DUAL DISTRIBUTION SUCCESSFUL\n");

			} else if (strcmp(part_type, "nodal") == 0) {
				part_key = 3;
				if (METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL,
						&nparts,
						NULL, NULL, objval, epart, npart) != METIS_OK)
					printf("ERROR IN METIS NODAL DISTRIBUTION\n");
				else
					printf("METIS NODAL DISTRIBUTION SUCCESSFUL\n");

			} else {
				part_key = 1;
				r = (nintcf_m + 1) % nprocs;
				q = (nintcf_m + 1) / nprocs;
				for (i = 0; i < nprocs; ++i) {
					if (i >= r)
						s = 0;
					for (j = 0; j < q + s; ++j) {
						//g_l[k] = j;
						epart[k++] = i;
					}
				}
			}

			for (i = 0; i <= nintcf_m; ++i) {
				nintcf_a[epart[i]]++;
				nextcf_a[epart[i]]++;
			}
			displs[0] = 0;
			displs2[0] = 0;
			for (i = 1; i < nprocs; ++i) {
				displs[i] = displs[i - 1] + nintcf_a[i - 1];
				displs2[i] = displs[i];
			}

			for (i = 0; i <= nintcf_m; ++i) {

				s = displs2[epart[i]]++;
				bp_m[s] = bp_m0[i];
				bs_m[s] = bs_m0[i];
				be_m[s] = be_m0[i];
				bn_m[s] = bn_m0[i];
				bw_m[s] = bw_m0[i];
				bl_m[s] = bl_m0[i];
				bh_m[s] = bh_m0[i];
				su_m[s] = su_m0[i];
				l_g[s] = i;
				g_l[i] = s - displs[epart[i]];
			}

			for (i = 0; i <= nextcf_m; ++i) {
				epart_ghost[i] = -1;
				g_l_ghost[i] = -1;
			}
			k = 0;
			for (i = 0; i <= nintcf_m; ++i) {
				for (j = 0; j < 6; ++j) {
					lcc_ij = lcc_m[i][j];
					if (lcc_ij > nintcf_m || epart[i] != epart[lcc_ij]) {
						if (epart_ghost[lcc_ij] != epart[i]) {
							epart_ghost[lcc_ij] = epart[i];
							g_l_ghost[lcc_ij] = ++nextcf_a[epart[i]];
						}
						lcc_a_m[k++] = g_l_ghost[lcc_ij];
					} else {
						lcc_a_m[k++] = g_l[lcc_ij];
					}
				}
			}
			nintcf_a2[0] = 6 * nintcf_a[0];
			displs2[0] = 0;
			for (i = 0; i < nprocs; ++i) {
				nintcf_a2[i] = 6 * nintcf_a[i];
				displs2[i] = displs2[i - 1] + nintcf_a2[i - 1];
			}

		}
		*elems = (int*) calloc(sizeof(int), 1);
		*points = (int**) calloc(sizeof(int*), 1);
		*points_count = 0;

		MPI_Scatter(nintcf_a, 1, MPI_INT, nintcf, 1, MPI_INT, 0,
		MPI_COMM_WORLD);
		MPI_Scatter(nextcf_a, 1, MPI_INT, nextcf, 1, MPI_INT, 0,
		MPI_COMM_WORLD);

		*var = (double*) calloc(sizeof(double), *nextcf);
		*cgup = (double*) calloc(sizeof(double), *nextcf);
		*cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

		*bp = (double*) malloc(sizeof(double) * (*nextcf));
		*bs = (double*) malloc(sizeof(double) * (*nextcf));
		*be = (double*) malloc(sizeof(double) * (*nextcf));
		*bn = (double*) malloc(sizeof(double) * (*nextcf));
		*bw = (double*) malloc(sizeof(double) * (*nextcf));
		*bl = (double*) malloc(sizeof(double) * (*nextcf));
		*bh = (double*) malloc(sizeof(double) * (*nextcf));
		*su = (double*) malloc(sizeof(double) * (*nextcf));

		*local_global_index = (int*) malloc(sizeof(int) * (*nintcf));
		lcc_a = (int*) malloc(sizeof(int) * 6 * (*nintcf));

		/*********************SCATTER ALL THE VALUES TO THE OTHER PROCESSES***********************/

		MPI_Scatterv(bp_m, nintcf_a, displs, MPI_DOUBLE, *bp, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Scatterv(bs_m, nintcf_a, displs, MPI_DOUBLE, *bs, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(be_m, nintcf_a, displs, MPI_DOUBLE, *be, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bn_m, nintcf_a, displs, MPI_DOUBLE, *bn, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bw_m, nintcf_a, displs, MPI_DOUBLE, *bw, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bl_m, nintcf_a, displs, MPI_DOUBLE, *bl, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bh_m, nintcf_a, displs, MPI_DOUBLE, *bh, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);

		MPI_Scatterv(su_m, nintcf_a, displs, MPI_DOUBLE, *su, *nintcf,
		MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(l_g, nintcf_a, displs, MPI_INT, *local_global_index,
				*nintcf, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(l_g, nintcf_a, displs, MPI_INT, *local_global_index,
				*nintcf, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(lcc_a_m, nintcf_a2, displs2, MPI_INT, lcc_a, 6 * (*nintcf),
		MPI_INT, 0, MPI_COMM_WORLD);

		*lcc = (int**) malloc((*nintcf) * sizeof(int*));
		for (i = 0; i < *nintcf; ++i) {
			(*lcc)[i] = (int *) malloc(6 * sizeof(int));
		}
		for (i = 0; i < *nintcf; ++i) {
			for (j = 0; j < 6; ++j) {
				(*lcc)[i][j] = lcc_a[6 * i + j];
			}
		}

		(*nextcf)--;
		(*nintcf)--;
		*nintci = 0;
		/**********************INITIALIZE THE VARIABLES*****************************************/
		for (i = 0; i <= 10; i++) {
			(*cnorm)[i] = 1.0;
		}
		for (i = (*nintci); i <= (*nintcf); i++) {
			(*var)[i] = 0.0;
		}
		for (i = (*nintci); i <= (*nintcf); i++) {
			(*cgup)[i] = 1.0 / ((*bp)[i]);
		}
		for (i = (*nintcf) + 1; i <= (*nextcf); i++) {
			(*var)[i] = 0.0;
			(*cgup)[i] = 0.0;
			(*bs)[i] = 0.0;
			(*be)[i] = 0.0;
			(*bn)[i] = 0.0;
			(*bw)[i] = 0.0;
			(*bh)[i] = 0.0;
			(*bl)[i] = 0.0;
		}

//		/*********** TEST DISTRIBUTION VTK OUTPUT FOR SU**************/
//		char file_out_su[100];
//		sprintf(file_out_su, "%s.SU.%s.rank%d.vtk", file_in, part_type, myrank);
//		test_distribution(file_in, file_out_su, *local_global_index,
//				*nintcf + 1, *su);
//
//		/*********** TEST DISTRIBUTION VTK OUTPUT FOR CGUP**************/
//		char file_out_cgup[100];
//		sprintf(file_out_cgup, "%s.CGUP.%s.rank%d.vtk", file_in, part_type,
//				myrank);
//		test_distribution(file_in, file_out_cgup, *local_global_index,
//				*nintcf + 1, *cgup);
//		/*********** TEST DISTRIBUTION FOR PARTITIONS**************************/
//		int ne_l = *nintcf + 1;
//		double myrank_a[ne_l];
//		for (i = 0; i < ne_l; i++) {
//			myrank_a[i] = myrank;
//		}
//		char file_out_part[100];
//		sprintf(file_out_part, "%s.%s.rank%d.vtk", file_in, part_type, myrank);
//		test_distribution(file_in, file_out_part, *local_global_index, ne_l,
//				myrank_a);

	} else if (strcmp(read_type, "allread") == 0) {
		read_key = 2;
//ALL READ APPROACH********************
		int nintcf_m, nextci_m, nextcf_m;

		double *bp_m0, *bs_m0, *be_m0, *bn_m0, *bw_m0, *bl_m0, *bh_m0, *su_m0;
		int r, q, k = 0;

		int f_status = read_binary_geo(file_in, nintci, &nintcf_m, &nextci_m,
				&nextcf_m, &lcc_m, &bs_m0, &be_m0, &bn_m0, &bw_m0, &bl_m0,
				&bh_m0, &bp_m0, &su_m0, points_count, points, elems);
		if (f_status != 0)
			return f_status;

		r = (nintcf_m + 1) % nprocs;
		q = (nintcf_m + 1) / nprocs;

		int temp1 = myrank * q + MIN(myrank, r);
		int temp2 = (myrank + 1) * q + MIN(myrank + 1, r) - 1;

		*nintci = 0;
		*nintcf = temp2 - temp1;

		*bp = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*su = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*bs = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*be = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*bn = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*bw = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*bl = (double*) malloc(sizeof(double) * (*nintcf + 1));
		*bh = (double*) malloc(sizeof(double) * (*nintcf + 1));

		idx_t ne = nintcf_m + 1;
		idx_t nn = *points_count;
		idx_t ncommon = 4;
		idx_t nparts = nprocs;

		idx_t *eptr = (idx_t*) malloc((ne + 1) * sizeof(idx_t));
		idx_t *eind = (idx_t*) malloc(ne * 8 * sizeof(idx_t));
		idx_t *objval = (idx_t*) malloc(1 * sizeof(idx_t));
		idx_t *epart = (idx_t*) malloc((ne) * sizeof(idx_t));
		idx_t *npart = (idx_t*) malloc((nn) * sizeof(idx_t));

		for (i = 0; i < (ne + 1); i++) {
			eptr[i] = i * 8;
		}

		for (i = 0; i < (ne * 8); i++) {
			eind[i] = (*elems)[i];
		}

		/**********METIS DISTRIBUTION*************/
		if (strcmp(part_type, "dual") == 0) {
			part_key = 2;
			if (METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon,
					&nparts, NULL, NULL, objval, epart, npart) != METIS_OK)
				printf("ERROR IN METIS DUAL DISTRIBUTION rank %d\n", myrank);
			else
				printf("METIS DUAL DISTRIBUTION SUCCESSFUL\n");

			int counter = 0;

			for (i = 0; i < ne; i++) {
				counter += (epart[i] == myrank);
			}

			*local_global_index = (int*) malloc(sizeof(int) * (counter + 1));

			k = 0;
			for (i = 0; i < ne; i++) {
				if (epart[i] == myrank) {

					(*local_global_index)[k] = i;
					(*bp)[k] = bp_m0[i];
					(*bs)[k] = bp_m0[i];
					(*be)[k] = bp_m0[i];
					(*bn)[k] = bp_m0[i];
					(*bw)[k] = bp_m0[i];
					(*bl)[k] = bp_m0[i];
					(*bh)[k] = bp_m0[i];
					(*su)[k] = su_m0[i];
					k++;
				}
			}

		} else if (strcmp(part_type, "nodal") == 0) {
			part_key = 3;
			if (METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nparts,
			NULL, NULL, objval, epart, npart) != METIS_OK)
				printf("ERROR IN METIS NODAL DISTRIBUTION\n");
			else
				printf("METIS NODAL DISTRIBUTION SUCCESSFUL\n");

			int counter = 0;

			for (i = 0; i < ne; i++) {
				counter += (epart[i] == myrank);
			}
			*local_global_index = (int*) malloc(sizeof(int) * (counter));

			k = 0;
			for (i = 0; i < ne; i++) {
				if (epart[i] == myrank) {
					(*local_global_index)[k] = i;
					(*bp)[k] = bp_m0[i];
					(*bs)[k] = bp_m0[i];
					(*be)[k] = bp_m0[i];
					(*bn)[k] = bp_m0[i];
					(*bw)[k] = bp_m0[i];
					(*bl)[k] = bp_m0[i];
					(*bh)[k] = bp_m0[i];
					(*su)[k] = su_m0[i];
					k++;
				}
			}

		} else {
			part_key = 1;
			*local_global_index = (int*) malloc(sizeof(int) * (*nintcf + 1));
			k = 0;
			for (i = *nintci; i <= *nintcf; i++) {
				(*local_global_index)[k++] = i;
			}

			memcpy(*bp, &bp_m0[temp1], (*nintcf));
			memcpy(*bs, &bs_m0[temp1], (*nintcf));
			memcpy(*be, &be_m0[temp1], (*nintcf));
			memcpy(*bn, &bn_m0[temp1], (*nintcf));
			memcpy(*bw, &bw_m0[temp1], (*nintcf));
			memcpy(*bl, &bl_m0[temp1], (*nintcf));
			memcpy(*bh, &bh_m0[temp1], (*nintcf));
			memcpy(*su, &su_m0[temp1], (*nintcf));

		}
		*elems = (int*) calloc(sizeof(int), 1);
		*points = (int**) calloc(sizeof(int*), 1);
		*points_count = 0;

		*var = (double*) calloc(sizeof(double), (*nintcf + 1));
		*cgup = (double*) calloc(sizeof(double), (*nintcf + 1));
		*cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

		for (i = 0; i <= 10; i++) {
			(*cnorm)[i] = 1.0;
		}

		for (i = (*nintci); i <= (*nintcf); i++) {
			(*var)[i] = 0.0;
		}

		for (i = (*nintci); i <= (*nintcf); i++) {
			(*cgup)[i] = 1.0 / ((*bp)[i]);
		}

		for (i = (*nintcf) + 1; i <= (*nintcf); i++) {
			(*var)[i] = 0.0;
			(*cgup)[i] = 0.0;
			(*bs)[i] = 0.0;
			(*be)[i] = 0.0;
			(*bn)[i] = 0.0;
			(*bw)[i] = 0.0;
			(*bh)[i] = 0.0;
			(*bl)[i] = 0.0;
		}

//		/***********TEST DISTRIBUTION VTK OUTPUT FOR SU**************/
//		char file_out_su[100];
//		sprintf(file_out_su, "%s.SU.%s.rank%d.vtk", file_in, part_type, myrank);
//		test_distribution(file_in, file_out_su, *local_global_index,
//				*nintcf + 1, *su);
//
//		/***********TEST DISTRIBUTION VTK OUTPUT FOR CGUP**************/
//		char file_out_cgup[100];
//		sprintf(file_out_cgup, "%s.CGUP.%s.rank%d.vtk", file_in, part_type,
//				myrank);
//		test_distribution(file_in, file_out_cgup, *local_global_index,
//				*nintcf + 1, *cgup);
//
//		/***********TEST DISTRIBUTION FOR PARTITIONS VTK**************************/
//		int ne_l = *nintcf + 1;
//		double myrank_a[ne_l];
//		for (i = 0; i < ne_l; i++) {
//			myrank_a[i] = myrank;
//		}
//		char file_out_part[100];
//		sprintf(file_out_part, "%s.%s.rank%d.vtk", file_in, part_type, myrank);
//		test_distribution(file_in, file_out_part, *local_global_index,
//				ne_l, myrank_a);

	}
	t2 = PAPI_get_virt_usec();
	t = t2 - t1;
	write_pstats_exectime(input_key, part_key, read_key, myrank, t);

	local_intc = *nintcf + 1;
	local_extc = *nextcf - *nintcf;
	write_pstats_partition(input_key, part_key, myrank, local_intc, local_extc);
	return 0;
}
