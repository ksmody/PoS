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
		double** cnorm, int** local_global_index, int** global_local_index,
		int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst,
		int **recv_cnt, int*** recv_lst) {

//	if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT) {
//		printf("Error: PAPI_library_init!");
//		exit(1);
//	}
//	long long int t1, t2;
//	t1 = PAPI_get_virt_usec();
//	int read_key = 1;
	int input_key = 1, part_key = 1;
	int *nintcf_a, *nextcf_a, *displs, *nintcf_a2, *displs2, *nsend_a;
	int *epart_recv, *epart_recv_m, *epart_send, *epart_send_m, *sendl, *sendl_m;
	double *bs_m, *be_m, *bn_m, *bw_m, *bl_m, *bh_m, *bp_m, *su_m;
	double *bp_m0, *bs_m0, *be_m0, *bn_m0, *bw_m0, *bl_m0, *bh_m0, *su_m0;
	int *lcc_a_m, *lcc_a, *l_g, *g_l, *g_l_2, *epart2;
	int *rank_ngh, *send_pos, *recv_pos;
	int nintcf_m, nextci_m, nextcf_m, ngh, nextc, nsend;
	int **lcc_m, **points_m, *elems_m;
	int i, j, k = 0, n, p, q, r, s = 1, t, i2 = 0, u[6];
	idx_t ne, nn, ncommon, nparts, *eptr, *eind, *objval, *npart, *epart;

	if (strcmp(file_in, "tjunc") == 0)
		input_key = 1;
	else if (strcmp(file_in, "drall") == 0)
		input_key = 2;
	else if (strcmp(file_in, "pent") == 0)
		input_key = 3;
	else if (strcmp(file_in, "cojack") == 0)
		input_key = 4;

	if ((strcmp(read_type, "oneread") == 0 && myrank == 0)
			|| strcmp(read_type, "allread") == 0) {

		int f_status = read_binary_geo(file_in, nintci, &nintcf_m, &nextci_m,
				&nextcf_m, &lcc_m, &bs_m0, &be_m0, &bn_m0, &bw_m0, &bl_m0,
				&bh_m0, &bp_m0, &su_m0, points_count, &points_m, &elems_m);
		if (f_status != 0)
			return f_status;

		g_l = (int*) malloc(sizeof(int) * (nintcf_m + 1));
		g_l_2 = (int*) malloc(sizeof(int) * (nextcf_m + 1));
		epart2 = (int*) malloc(sizeof(int) * (nextcf_m + 1));

		lcc_a_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
		sendl_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
		epart_recv_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
		epart_send_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
		l_g = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);


		bp_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		bs_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		be_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		bn_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		bw_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		bl_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		bh_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
		su_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));

		displs = (int*) calloc(sizeof(int), nprocs);
		displs2 = (int*) calloc(sizeof(int), nprocs);
		nintcf_a = (int*) calloc(sizeof(int), nprocs);
		nintcf_a2 = (int*) calloc(sizeof(int), nprocs);
		nextcf_a = (int*) calloc(sizeof(int), nprocs);
		nsend_a = (int*) calloc(sizeof(int), nprocs);


		ne = nintcf_m + 1;
		nn = *points_count;
		ncommon = 4;
		nparts = nprocs;

		eptr = (idx_t*) malloc((ne + 1) * sizeof(idx_t*));
		eind = (idx_t*) malloc(ne * 8 * sizeof(idx_t*));
		objval = (idx_t*) malloc(1 * sizeof(idx_t*));
		npart = (idx_t*) malloc((nn) * sizeof(idx_t*));
		epart = (idx_t*) malloc((ne) * sizeof(idx_t*));


		for (i = 0; i < (ne + 1); i++) {
			eptr[i] = i * 8;
		}
		for (i = 0; i < (ne * 8); i++) {
			eind[i] = elems_m[i];
		}

		/**********METIS DISTRIBUTION*************/
		if (strcmp(part_type, "dual") == 0) {
			part_key = 2;
			if (METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon,
					&nparts, NULL, NULL, objval, epart, npart) != METIS_OK)
				printf("ERROR IN METIS DUAL DISTRIBUTION\n");
			else
				printf("METIS DUAL DISTRIBUTION SUCCESSFUL\n");

		} else if (strcmp(part_type, "nodal") == 0) {
			part_key = 3;
			if (METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nparts,
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
					epart[k++] = i;
				}
			}
			printf("CLASSIC DISTRIBUTION SUCCESSFUL\n");
		}

		for (i = 0; i <= nintcf_m; ++i) {
			nintcf_a[epart[i]]++;
			nextcf_a[epart[i]]++;
		}

		for (i = 1; i < nprocs; ++i) {
			displs[i] = displs[i - 1] + nintcf_a[i - 1];
			displs2[i] = displs[i];
		}

		for (i = 0; i <= nintcf_m; ++i) {
			s = displs[epart[i]];
			k = displs2[epart[i]]++;
			bp_m[k] = bp_m0[i];
			bs_m[k] = bs_m0[i];
			be_m[k] = be_m0[i];
			bn_m[k] = bn_m0[i];
			bw_m[k] = bw_m0[i];
			bl_m[k] = bl_m0[i];
			bh_m[k] = bh_m0[i];
			su_m[k] = su_m0[i];
			g_l[i] = k - s;
			l_g[5 * s + k] = i;
		}

		for (i = 0; i < nprocs; ++i) {
			nintcf_a2[i] = 6 * nintcf_a[i];
			displs2[i] = 6 * displs[i];
		}

		for (i = 0; i <= nextcf_m; ++i) {
			epart2[i] = -1;
			g_l_2[i] = -1;
		}
		k = 0;
		for (p = 0; p < nprocs; ++p) {
			s = displs2[p];
			n = nextcf_a[p];
			for (i = 0; i < nintcf_a[p]; ++i) {
				i2 = l_g[s + i];
				for (j = 0; j < 6; ++j)
					u[j] = -1;
				for (j = 0; j < 6; ++j) {
					ngh = lcc_m[i2][j];
					if (ngh > nintcf_m || p != epart[ngh]) {
						if (epart2[ngh] != p) {
							epart2[ngh] = p;
							g_l_2[ngh] = n;
							l_g[s + n] = ngh;
							epart_recv_m[s + n - nintcf_a[p]] =
									ngh > nintcf_m ? -1 : epart[ngh];
							nextcf_a[p]++;
							n++;
						}
						if (ngh <= nintcf_m) {
							t = epart[ngh];
							if (t != u[0] && t != u[1] && t != u[2] && t != u[3]
									&& t != u[4] && t != u[5]) {
								u[j] = t;
								sendl_m[s + nsend_a[p]] = i;
								epart_send_m[s + nsend_a[p]] = t;
								nsend_a[p]++;
							}
						}
						lcc_a_m[k++] = g_l_2[ngh];
					} else {
						lcc_a_m[k++] = g_l[ngh];
					}
				}
			}
		}
	}
	if (strcmp(read_type, "oneread") == 0) {
		//read_key = 1;
		MPI_Scatter(nintcf_a, 1, MPI_INT, nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(nextcf_a, 1, MPI_INT, nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(nsend_a, 1, MPI_INT, &nsend, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		//read_key = 2;
		*nintcf = nintcf_a[myrank];
		*nextcf = nextcf_a[myrank];
		nsend = nsend_a[myrank];
	}
	nextc = *nextcf - *nintcf;

	*points_count = 1;
	*points = (int **) calloc(*points_count, sizeof(int *));
	for (i = 0; i < *points_count; ++i) {
		(*points)[i] = (int *) calloc(3, sizeof(int));
	}
	*elems = (int *) malloc(sizeof(int) * (*nintcf) * 8);

	*var = (double*) calloc(sizeof(double), *nextcf);
	*cgup = (double*) calloc(sizeof(double), *nextcf);
	*cnorm = (double*) calloc(sizeof(double), *nintcf);

	*bp = (double*) malloc(sizeof(double) * (*nextcf));
	*bs = (double*) malloc(sizeof(double) * (*nextcf));
	*be = (double*) malloc(sizeof(double) * (*nextcf));
	*bn = (double*) malloc(sizeof(double) * (*nextcf));
	*bw = (double*) malloc(sizeof(double) * (*nextcf));
	*bl = (double*) malloc(sizeof(double) * (*nextcf));
	*bh = (double*) malloc(sizeof(double) * (*nextcf));
	*su = (double*) malloc(sizeof(double) * (*nextcf));

	*local_global_index = (int*) malloc(sizeof(int) * (*nextcf));
	lcc_a = (int*) malloc(sizeof(int) * 6 * (*nintcf));
	sendl = (int*) malloc(sizeof(int) * nsend);
	epart_send = (int*) malloc(sizeof(int) * nsend);
	epart_recv = (int*) malloc(sizeof(int) * nextc);


	/*********************SCATTER ALL THE VALUES TO THE OTHER PROCESSES***********************/

	if (strcmp(read_type, "oneread") == 0) {
		MPI_Scatterv(bp_m, nintcf_a, displs, MPI_DOUBLE, *bp, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bs_m, nintcf_a, displs, MPI_DOUBLE, *bs, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(be_m, nintcf_a, displs, MPI_DOUBLE, *be, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bn_m, nintcf_a, displs, MPI_DOUBLE, *bn, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bw_m, nintcf_a, displs, MPI_DOUBLE, *bw, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bl_m, nintcf_a, displs, MPI_DOUBLE, *bl, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(bh_m, nintcf_a, displs, MPI_DOUBLE, *bh, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(su_m, nintcf_a, displs, MPI_DOUBLE, *su, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(l_g, nextcf_a, displs2, MPI_INT, *local_global_index, *nextcf,	MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(lcc_a_m, nintcf_a2, displs2, MPI_INT, lcc_a, 6 * (*nintcf), MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		for (i = 0; i < *nintcf; i++) {
			(*bp)[i] = bp_m[displs[myrank] + i];
			(*bn)[i] = bn_m[displs[myrank] + i];
			(*be)[i] = be_m[displs[myrank] + i];
			(*bw)[i] = bw_m[displs[myrank] + i];
			(*bs)[i] = bs_m[displs[myrank] + i];
			(*bl)[i] = bl_m[displs[myrank] + i];
			(*bh)[i] = bh_m[displs[myrank] + i];
			(*su)[i] = su_m[displs[myrank] + i];
		}
		k = 0;
		for (i = displs2[myrank]; i < displs2[myrank] + *nextcf; ++i) {
			(*local_global_index)[k++] = l_g[i];
		}
		k = 0;
		for (i = displs2[myrank]; i < displs2[myrank] + 6 * (*nintcf); ++i) {
			lcc_a[k++] = lcc_a_m[i];
		}
	}

	*lcc = (int**) malloc((*nintcf) * sizeof(int*));
	for (i = 0; i < *nintcf; ++i) {
		(*lcc)[i] = (int *) malloc(6 * sizeof(int));
	}
	for (i = 0; i < *nintcf; ++i) {
		for (j = 0; j < 6; ++j) {
			(*lcc)[i][j] = lcc_a[6 * i + j];
		}
	}

	if (strcmp(read_type, "oneread") == 0) {
		if (myrank == 0) {
			for (i = 0; i < nprocs; ++i) {
				nextcf_a[i] = nextcf_a[i] - nintcf_a[i];
			}
		}
		MPI_Scatterv(sendl_m, nsend_a, displs2, MPI_INT, sendl, nsend, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(epart_send_m, nsend_a, displs2, MPI_INT, epart_send, nsend, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(epart_recv_m, nextcf_a, displs2, MPI_INT, epart_recv, nextc, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		for (i = 0; i < nsend; ++i) {
			sendl[i] = sendl_m[displs2[myrank] + i];
			epart_send[i] = epart_send_m[displs2[myrank] + i];
		}
		for (i = 0; i < nextc; ++i) {
			epart_recv[i] = epart_recv_m[displs2[myrank] + i];
		}
	}

	rank_ngh = (int *) malloc(sizeof(int) * nprocs);
	for (i = 0; i < nprocs; ++i) {
		rank_ngh[i] = -1;
	}
	*nghb_cnt = 0;
	for (i = 0; i < nextc; ++i) {
		j = epart_recv[i];
		if (j != -1) {
			if (rank_ngh[j] == -1) {
				rank_ngh[j] = (*nghb_cnt)++;
			}
		}
	}
	send_pos = (int*) calloc(sizeof(int), *nghb_cnt);
	recv_pos = (int*) calloc(sizeof(int), *nghb_cnt);
	*send_cnt = (int*) calloc(sizeof(int), *nghb_cnt);
	*recv_cnt = (int*) calloc(sizeof(int), *nghb_cnt);
	*send_lst = (int**) calloc(sizeof(int*), *nghb_cnt);
	*recv_lst = (int**) calloc(sizeof(int*), *nghb_cnt);
	*nghb_to_rank = (int*) calloc(sizeof(int), *nghb_cnt);


	for (i = 0; i < nsend; ++i) {
		epart_send[i] = rank_ngh[epart_send[i]];
		(*send_cnt)[epart_send[i]]++;
	}
	for (i = 0; i < nextc; ++i) {
		j = epart_recv[i];
		if (j != -1) {
			k = rank_ngh[j];
			(*nghb_to_rank)[k] = j;
			(*recv_cnt)[k]++;
			epart_recv[i] = k;
		}
	}
	for (i = 0; i < *nghb_cnt; ++i) {
		(*send_lst)[i] = (int *) malloc(sizeof(int) * (*send_cnt)[i]);
		(*recv_lst)[i] = (int *) malloc(sizeof(int) * (*recv_cnt)[i]);
	}

	for (i = 0; i < nsend; ++i) {
		j = epart_send[i];
		(*send_lst)[j][send_pos[j]++] = sendl[i];
	}
	for (i = 0; i < nextc; ++i) {
		j = epart_recv[i];
		if (j != -1) {
			(*recv_lst)[j][recv_pos[j]++] = *nintcf + i;
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

	/*********** TEST DISTRIBUTION FOR SEND & RECEIVE LISTS**************************/

	n = *nintcf + 1;
	int l_g_no_boundary[*nextcf + 1];
	double myrank_s[*nextcf + 1];
	double myrank_r[*nextcf + 1];
	for (i = 0; i < n; ++i) {
		l_g_no_boundary[i] = (*local_global_index)[i];
		myrank_s[i] = myrank + 1;
		myrank_r[i] = myrank + 1;
	}
	for (i = 0; i < *nghb_cnt; ++i) {
		for (j = 0; j < (*send_cnt)[i]; ++j) {
			myrank_s[(*send_lst)[i][j]] = (*nghb_to_rank)[i] + 1;
		}
	}
	for (i = 0; i < *nghb_cnt; ++i) {
		for (j = 0; j < (*recv_cnt)[i]; ++j) {
			l_g_no_boundary[n] = (*local_global_index)[(*recv_lst)[i][j]];
			myrank_r[n] = (*nghb_to_rank)[i] + 1;
			n++;
		}
	}
	char file_out_send[100];
	sprintf(file_out_send, "%s.%s.SEND.rank%d.vtk", file_in, part_type, myrank);
	test_distribution(file_in, file_out_send, l_g_no_boundary, *nintcf + 1, myrank_s);
	char file_out_recv[100];
	sprintf(file_out_recv, "%s.%s.RECV.rank%d.vtk", file_in, part_type, myrank);
	test_distribution(file_in, file_out_recv, l_g_no_boundary, n, myrank_r);

	for (i = 0; i < *nghb_cnt; ++i) {
		write_pstats_communication(input_key, part_key, myrank, nprocs,
				*nghb_cnt, i, *send_cnt, *send_lst, *recv_cnt, *recv_lst);
	}
//	t2 = PAPI_get_virt_usec();
//	write_pstats_exectime(input_key, part_key, read_key, myrank, t2 - t1);
//	write_pstats_partition(input_key, part_key, myrank, *nintcf + 1, nextc);

	if ((strcmp(read_type, "oneread") == 0 && myrank == 0)
			|| strcmp(read_type, "allread") == 0) {
		free(g_l);
		free(g_l_2);
		free(epart2);
		free(lcc_a_m);
		free(sendl_m);
		free(epart_recv_m);
		free(epart_send_m);
		free(l_g);
		free(bp_m);
		free(bs_m);
		free(be_m);
		free(bn_m);
		free(bw_m);
		free(bl_m);
		free(bh_m);
		free(su_m);
		free(displs);
		free(displs2);
		free(nintcf_a);
		free(nintcf_a2);
		free(nextcf_a);
		free(nsend_a);
		free(eptr);
		free(eind);
		free(objval);
		free(npart);
		free(epart);
	}
	free(lcc_a);
	free(sendl);
	free(epart_send);
	free(epart_recv);
	free(rank_ngh);
	free(send_pos);
	free(recv_pos);

	return 0;
}

