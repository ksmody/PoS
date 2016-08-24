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

void partition(char *part_type, idx_t nprocs, idx_t nn, idx_t ne, int *elems, idx_t **epart) {
	idx_t ncommon = 4;
	idx_t *eptr = (idx_t*) malloc((ne + 1) * sizeof(idx_t*));
	idx_t *eind = (idx_t*) malloc(ne * 8 * sizeof(idx_t*));
	idx_t *objval = (idx_t*) malloc(1 * sizeof(idx_t*));
	idx_t *npart = (idx_t*) malloc((nn) * sizeof(idx_t*));
	*epart = (idx_t*) malloc((ne) * sizeof(idx_t*));
	/* epart[i] - rank of the i-th element */
	int i, j, k = 0, s = 1;

	/* Metis parameters - elements and their displacements */
	for (i = 0; i < (ne * 8); i++) {
		eind[i] = elems[i];
	}
	for (i = 0; i < (ne + 1); i++) {
		eptr[i] = i * 8;
	}
	/* Domain distribution - epart definition, depending on the partition type */
	if (strcmp(part_type, "dual") == 0) {
		if (METIS_PartMeshDual(&ne, &nn, eptr, eind, NULL, NULL, &ncommon,
				&nprocs, NULL, NULL, objval, *epart, npart) != METIS_OK)
			printf("ERROR IN METIS DUAL DISTRIBUTION\n");
		else
			printf("METIS DUAL DISTRIBUTION SUCCESSFUL\n");

	} else if (strcmp(part_type, "nodal") == 0) {
		if (METIS_PartMeshNodal(&ne, &nn, eptr, eind, NULL, NULL, &nprocs,
				NULL, NULL, objval, *epart, npart) != METIS_OK)
			printf("ERROR IN METIS NODAL DISTRIBUTION\n");
		else
			printf("METIS NODAL DISTRIBUTION SUCCESSFUL\n");
	} else if (strcmp(part_type, "classic") == 0){
		/* i-th chunk of the elements goes to the rank i */
		for (i = 0; i < nprocs; ++i) {
			if (i >= ne % nprocs)
				s = 0;
			for (j = 0; j < ne / nprocs + s; ++j) {
				(*epart)[k++] = i;
			}
		}
		printf("CLASSIC DISTRIBUTION SUCCESSFUL\n");
	}
	free(eptr);
	free(eind);
	free(objval);
	free(npart);
}

void arrangeData(int nintcf_m, int **nintcf_a, int **nextcf_a, int **nintcf_a2, int **displs, int **displs2,
		double **bs_m0, double **be_m0, double **bn_m0, double **bw_m0, double **bl_m0, double **bh_m0,
		double **bs_m, double **be_m, double **bn_m, double **bw_m, double **bl_m, double **bh_m,
		double **bp_m0, double **su_m0, double **bp_m, double **su_m,
		int **l_g, int **g_l, idx_t *epart, int nprocs) {
	int i, s, k;
	*g_l = (int*) malloc(sizeof(int) * (nintcf_m + 1));
	*l_g = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 7);

	*bp_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*bs_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*be_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*bn_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*bw_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*bl_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*bh_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));
	*su_m = (double*) malloc(sizeof(double) * (nintcf_m + 1));

	*displs = (int*) calloc(sizeof(int), nprocs);
	*displs2 = (int*) calloc(sizeof(int), nprocs);
	*nintcf_a = (int*) calloc(sizeof(int), nprocs);
	*nintcf_a2 = (int*) calloc(sizeof(int), nprocs);
	*nextcf_a = (int*) calloc(sizeof(int), nprocs);

	/* Number of (internal) elements for each process */
	for (i = 0; i <= nintcf_m; ++i) {
		(*nintcf_a)[epart[i]]++;
		(*nextcf_a)[epart[i]]++;
	}
	/* Displacement of elements for each process */
	for (i = 1; i < nprocs; ++i) {
		(*displs)[i] = (*displs)[i - 1] + (*nintcf_a)[i - 1];
		(*displs2)[i] = (*displs)[i];
	}
	/* Rearranging of data before its distribution to other processes
	 * Definition of mappings between global and local indices for internal cells */
	for (i = 0; i <= nintcf_m; ++i) {
		s = (*displs)[epart[i]];
		k = (*displs2)[epart[i]]++;
		(*bp_m)[k] = (*bp_m0)[i];
		(*bs_m)[k] = (*bs_m0)[i];
		(*be_m)[k] = (*be_m0)[i];
		(*bn_m)[k] = (*bn_m0)[i];
		(*bw_m)[k] = (*bw_m0)[i];
		(*bl_m)[k] = (*bl_m0)[i];
		(*bh_m)[k] = (*bh_m0)[i];
		(*su_m)[k] = (*su_m0)[i];
		(*g_l)[i] = k - s;
		(*l_g)[6 * s + k] = i;
	}
	/* 6x number of elements and their displacement for each process */
	for (i = 0; i < nprocs; ++i) {
		(*nintcf_a2)[i] = 6 * (*nintcf_a)[i];
		(*displs2)[i] = 6 * (*displs)[i];
	}
	free(*bp_m0);
	free(*bs_m0);
	free(*be_m0);
	free(*bn_m0);
	free(*bw_m0);
	free(*bl_m0);
	free(*bh_m0);
	free(*su_m0);
}

void defineGhost(int nprocs, idx_t **epart, int *displs, int *displs2, int nintcf_m, int nextcf_m,
		int **epart_send_m, int **epart_recv_m, int **sendl_m, int **nsend_a, int *nintcf_a, int *nextcf_a,
		int **g_l, int *l_g, int ***lcc_m, int **lcc_a_m) {

	int i, i2, j, k = 0, p, q, r, s, t, n, ngh;
	/* Temporary arrays for checking the ghost cell replication */
	int *g_l2 = (int*) malloc(sizeof(int) * (nextcf_m + 1));
	int *epart2 = (int*) malloc(sizeof(int) * (nextcf_m + 1));
	for (i = 0; i <= nextcf_m; ++i) {
		epart2[i] = -1;
		g_l2[i] = -1;
	}
	*lcc_a_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
	*sendl_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
	*epart_recv_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
	*epart_send_m = (int*) malloc(sizeof(int) * (nintcf_m + 1) * 6);
	*nsend_a = (int*) calloc(sizeof(int), nprocs);
	/* Definition of ghost and real boundary cells for all of the partitions
	 * Redefinition of lcc array and its rearrangement before scattering
	 * Definition of local_global array for external cells
	 * Send_list creation */
	for (p = 0; p < nprocs; ++p) {
		q = displs[p];
		s = displs2[p];
		r = nintcf_a[p];
		n = nextcf_a[p];
		for (i = 0; i < r; ++i) {
			i2 = l_g[s + q + i];
			for (j = 0; j < 6; ++j) {
				ngh = (*lcc_m)[i2][j];
				/* Define an external cell as either a real boundary or a ghost cell */
				if (ngh > nintcf_m || p != (*epart)[ngh]) {
					/* Check if a ghost cell has been already defined */
					if (epart2[ngh] != p) {
						epart2[ngh] = p;
						g_l2[ngh] = n;
						l_g[s + q + n] = ngh;
						/* -1 distingushes between real boundaries and ghost cells */
						(*epart_recv_m)[s + n - r] = -1;
						/* Define a sendlist and its corresponding epart array, as well as
						 * the one for recvlist (ghost cells) for internal cells only */
						if (ngh <= nintcf_m) {
							t = (*epart)[ngh];
							(*epart_recv_m)[s + n - r] = t;
							(*epart_send_m)[displs2[t] + (*nsend_a)[t]] = p;
							(*sendl_m)[displs2[t] + (*nsend_a)[t]] = (*g_l)[ngh];
							(*nsend_a)[t]++;
						}
						nextcf_a[p]++; n++;
					}
					(*lcc_a_m)[k++] = g_l2[ngh];
				} else {
					(*lcc_a_m)[k++] = (*g_l)[ngh];
				}
			}
		}
	}
    for ( i = 0; i < nintcf_m + 1; i++ ) {
        free((*lcc_m)[i]);
    }
    free(*lcc_m);
	free(*epart);
	free(epart2);
	free(*g_l);
	free(g_l2);
}

void allocateSpace(char *read_type, int myrank, int *nintcf_a, int *nextcf_a, int* nintcf, int* nextcf,
		double** bs, double** be, double** bn, double** bw, double** bl, double** bh, double** bp,
		double** su, double** var, double** cgup, double** cnorm, int*** lcc, int **lcc_a,
		int* points_count, int*** points, int** elems, int** local_global_index,
		int *nsend_a, int *nsend, int **sendl, int **epart_send, int **epart_recv) {

	/* Distribution of the size of internal, external and sendlist elements */
	if (strcmp(read_type, "oneread") == 0) {
		MPI_Scatter(nintcf_a, 1, MPI_INT, nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(nextcf_a, 1, MPI_INT, nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatter(nsend_a, 1, MPI_INT, nsend, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		*nintcf = nintcf_a[myrank];
		*nextcf = nextcf_a[myrank];
		*nsend = nsend_a[myrank];
	}
	int i, nextc = *nextcf - *nintcf;

	if (strcmp(read_type, "oneread") == 0 && myrank != 0) {
		*points_count = 1;
		*points = (int **) calloc(*points_count, sizeof(int *));
		for (i = 0; i < *points_count; ++i) {
			(*points)[i] = (int *) calloc(3, sizeof(int));
		}
		*elems = (int *) malloc(sizeof(int) * 8);
	}

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
	*lcc_a = (int*) malloc(sizeof(int) * 6 * (*nintcf));
	*sendl = (int*) malloc(sizeof(int) * *nsend);
	*epart_send = (int*) malloc(sizeof(int) * *nsend);
	*epart_recv = (int*) malloc(sizeof(int) * nextc);
}

void scatterData(char *read_type, int myrank, int **nintcf_a, int **nextcf_a, int* nintcf, int* nextcf,
		double** bs, double** be, double** bn, double** bw, double** bl, double** bh, double** bp,
		double** bs_m, double** be_m, double** bn_m, double** bw_m, double** bl_m, double** bh_m, double** bp_m,
		double** su, int *lcc_a, double** su_m, int **lcc_a_m, int** local_global_index,
		int **nsend_a, int nsend, int *sendl, int *epart_send, int *epart_recv,
		int **sendl_m, int **epart_send_m, int **epart_recv_m, int **l_g,
		int **displs, int **displs2, int **nintcf_a2, int nprocs) {

	int i, k, nextc = *nextcf - *nintcf;

	/* Distribution of all the arrays to the remaining processes */
	if (strcmp(read_type, "oneread") == 0) {
		MPI_Scatterv(*bp_m, *nintcf_a, *displs, MPI_DOUBLE, *bp, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*bs_m, *nintcf_a, *displs, MPI_DOUBLE, *bs, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*be_m, *nintcf_a, *displs, MPI_DOUBLE, *be, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*bn_m, *nintcf_a, *displs, MPI_DOUBLE, *bn, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*bw_m, *nintcf_a, *displs, MPI_DOUBLE, *bw, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*bl_m, *nintcf_a, *displs, MPI_DOUBLE, *bl, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*bh_m, *nintcf_a, *displs, MPI_DOUBLE, *bh, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*su_m, *nintcf_a, *displs, MPI_DOUBLE, *su, *nintcf, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*lcc_a_m, *nintcf_a2, *displs2, MPI_INT, lcc_a, 6 * (*nintcf), MPI_INT, 0, MPI_COMM_WORLD);
		if (myrank == 0) {
			for (i = 0; i < nprocs; ++i) {
				(*nintcf_a2)[i] = 7 * (*displs)[i];
			}
		}
		MPI_Scatterv(*l_g, *nextcf_a, *nintcf_a2, MPI_INT, *local_global_index, *nextcf, MPI_INT, 0, MPI_COMM_WORLD);
		if (myrank == 0) {
			for (i = 0; i < nprocs; ++i) {
				(*nextcf_a)[i] = (*nextcf_a)[i] - (*nintcf_a)[i];
			}
		}
		MPI_Scatterv(*sendl_m, *nsend_a, *displs2, MPI_INT, sendl, nsend, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*epart_send_m, *nsend_a, *displs2, MPI_INT, epart_send, nsend, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Scatterv(*epart_recv_m, *nextcf_a, *displs2, MPI_INT, epart_recv, nextc, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else {
		for (i = 0; i < *nintcf; i++) {
			(*bp)[i] = (*bp_m)[(*displs)[myrank] + i];
			(*bn)[i] = (*bn_m)[(*displs)[myrank] + i];
			(*be)[i] = (*be_m)[(*displs)[myrank] + i];
			(*bw)[i] = (*bw_m)[(*displs)[myrank] + i];
			(*bs)[i] = (*bs_m)[(*displs)[myrank] + i];
			(*bl)[i] = (*bl_m)[(*displs)[myrank] + i];
			(*bh)[i] = (*bh_m)[(*displs)[myrank] + i];
			(*su)[i] = (*su_m)[(*displs)[myrank] + i];
		}
		k = 0;
		for (i = 7 * (*displs)[myrank]; i < 7 * (*displs)[myrank] + *nextcf; ++i) {
			(*local_global_index)[k++] = (*l_g)[i];
		}
		k = 0;
		for (i = (*displs2)[myrank]; i < (*displs2)[myrank] + 6 * (*nintcf); ++i) {
			lcc_a[k++] = (*lcc_a_m)[i];
		}
		for (i = 0; i < nsend; ++i) {
			sendl[i] = (*sendl_m)[(*displs2)[myrank] + i];
			epart_send[i] = (*epart_send_m)[(*displs2)[myrank] + i];
		}
		for (i = 0; i < nextc; ++i) {
			epart_recv[i] = (*epart_recv_m)[(*displs2)[myrank] + i];
		}
	}
	if (myrank == 0	|| strcmp(read_type, "allread") == 0) {
		free(*lcc_a_m);
		free(*sendl_m);
		free(*epart_recv_m);
		free(*epart_send_m);
		free(*l_g);
		free(*bp_m);
		free(*bs_m);
		free(*be_m);
		free(*bn_m);
		free(*bw_m);
		free(*bl_m);
		free(*bh_m);
		free(*su_m);
		free(*displs);
		free(*displs2);
		free(*nintcf_a);
		free(*nintcf_a2);
		free(*nextcf_a);
		free(*nsend_a);
	}
}

void initVar(int *nintci, int *nintcf, int *nextci, int *nextcf,
		int ***lcc, int **lcc_a, double *cnorm, double *var, double *cgup,
		double *bp, double *bs, double *be, double *bn, double *bw, double *bh, double *bl) {
	int i, j;
	*lcc = (int**) malloc((*nintcf) * sizeof(int*));
	for (i = 0; i < *nintcf; ++i) {
		(*lcc)[i] = (int *) malloc(6 * sizeof(int));
	}
	for (i = 0; i < *nintcf; ++i) {
		for (j = 0; j < 6; ++j) {
			(*lcc)[i][j] = (*lcc_a)[6 * i + j];
		}
	}
	free(*lcc_a);

	(*nextcf)--;
	(*nintcf)--;
	*nintci = 0;
	*nextci = *nintcf + 1;

	for (i = 0; i <= 10; i++) {
		cnorm[i] = 1.0;
	}
	for (i = (*nintci); i <= (*nintcf); i++) {
		var[i] = 0.0;
	}
	for (i = (*nintci); i <= (*nintcf); i++) {
		cgup[i] = 1.0 / (bp[i]);
	}
	for (i = (*nintcf) + 1; i <= (*nextcf); i++) {
		var[i] = 0.0;
		cgup[i] = 0.0;
		bs[i] = 0.0;
		be[i] = 0.0;
		bn[i] = 0.0;
		bw[i] = 0.0;
		bh[i] = 0.0;
		bl[i] = 0.0;
	}
}

void createList(int nsend, int nextc, int **epart_send, int **epart_recv, int **sendl, int nprocs, int ne,
		int *nghb_cnt, int **nghb_to_rank, int **send_cnt, int **recv_cnt, int ***send_lst, int ***recv_lst) {
	int i, j, k;
	int *rank_ngh = (int *) malloc(sizeof(int) * nprocs);
	for (i = 0; i < nprocs; ++i) {
		rank_ngh[i] = -1;
	}
	*nghb_cnt = 0;
	for (i = 0; i < nsend; ++i) {
		j = (*epart_send)[i];
		if (rank_ngh[j] == -1) {
			rank_ngh[j] = (*nghb_cnt)++;
		}
	}
	int *send_pos = (int*) calloc(sizeof(int), *nghb_cnt);
	int *recv_pos = (int*) calloc(sizeof(int), *nghb_cnt);
	*send_cnt = (int*) calloc(sizeof(int), *nghb_cnt);
	*recv_cnt = (int*) calloc(sizeof(int), *nghb_cnt);
	*send_lst = (int**) calloc(sizeof(int*), *nghb_cnt);
	*recv_lst = (int**) calloc(sizeof(int*), *nghb_cnt);
	*nghb_to_rank = (int*) calloc(sizeof(int), *nghb_cnt);

	for (i = 0; i < nsend; ++i) {
		(*epart_send)[i] = rank_ngh[(*epart_send)[i]];
		(*send_cnt)[(*epart_send)[i]]++;
	}
	for (i = 0; i < nextc; ++i) {
		j = (*epart_recv)[i];
		if (j != -1) {
			k = rank_ngh[j];
			(*nghb_to_rank)[k] = j;
			(*recv_cnt)[k]++;
			(*epart_recv)[i] = k;
		}
	}
	for (i = 0; i < *nghb_cnt; ++i) {
		(*send_lst)[i] = (int *) malloc(sizeof(int) * (*send_cnt)[i]);
		(*recv_lst)[i] = (int *) malloc(sizeof(int) * (*recv_cnt)[i]);
	}

	for (i = 0; i < nsend; ++i) {
		j = (*epart_send)[i];
		(*send_lst)[j][send_pos[j]++] = (*sendl)[i];
	}
	for (i = 0; i < nextc; ++i) {
		j = (*epart_recv)[i];
		if (j != -1) {
			(*recv_lst)[j][recv_pos[j]++] = ne + i;
		}
	}
	free(*sendl);
	free(*epart_send);
	free(*epart_recv);
	free(rank_ngh);
	free(send_pos);
	free(recv_pos);
}

int initialization(char* file_in, char* part_type, char* read_type, int nprocs,
		int myrank, int* nintci, int* nintcf, int* nextci, int* nextcf,
		int*** lcc, double** bs, double** be, double** bn, double** bw,
		double** bl, double** bh, double** bp, double** su, int* points_count,
		int*** points, int** elems, double** var, double** cgup, double** oc,
		double** cnorm, int** local_global_index, int** global_local_index,
		int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst,
		int **recv_cnt, int*** recv_lst) {

	double *bs_m, *be_m, *bn_m, *bw_m, *bl_m, *bh_m, *bp_m, *su_m;
	double *bp_m0, *bs_m0, *be_m0, *bn_m0, *bw_m0, *bl_m0, *bh_m0, *su_m0;

	int **lcc_m, *lcc_a_m, *lcc_a, *l_g, *g_l;

	int nintcf_m, nextci_m, nextcf_m, nsend;
	int *nintcf_a, *nextcf_a, *displs, *nintcf_a2, *displs2, *nsend_a;

	int *epart_recv, *epart_recv_m, *epart_send, *epart_send_m, *sendl, *sendl_m;
	idx_t *epart;

	if (myrank == 0 || strcmp(read_type, "allread") == 0) {

		int f_status = read_binary_geo(file_in, nintci, &nintcf_m, &nextci_m,
									   &nextcf_m, &lcc_m, &bs_m0, &be_m0, &bn_m0, &bw_m0, &bl_m0,
									   &bh_m0, &bp_m0, &su_m0, points_count, points, elems);
		if (f_status != 0)
			return f_status;

		partition(part_type, nprocs, *points_count, nintcf_m + 1, *elems, &epart);

		arrangeData(nintcf_m, &nintcf_a, &nextcf_a, &nintcf_a2, &displs, &displs2,
					&bs_m0, &be_m0, &bn_m0, &bw_m0, &bl_m0, &bh_m0,
					&bs_m, &be_m, &bn_m, &bw_m, &bl_m, &bh_m,
					&bp_m0, &su_m0, &bp_m, &su_m,
					&l_g, &g_l, epart, nprocs);

		defineGhost(nprocs, &epart, displs, displs2, nintcf_m, nextcf_m,
					&epart_send_m, &epart_recv_m, &sendl_m, &nsend_a, nintcf_a, nextcf_a,
					&g_l, l_g, &lcc_m, &lcc_a_m);
	}
	allocateSpace(read_type, myrank, nintcf_a, nextcf_a, nintcf, nextcf, bs, be, bn, bw, bl, bh, bp,
				  su, var, cgup, cnorm, lcc, &lcc_a, points_count, points, elems, local_global_index,
				  nsend_a, &nsend, &sendl, &epart_send, &epart_recv);

	scatterData(read_type, myrank, &nintcf_a, &nextcf_a, nintcf, nextcf, bs, be, bn, bw, bl, bh, bp,
				&bs_m, &be_m, &bn_m, &bw_m, &bl_m, &bh_m, &bp_m, su, lcc_a, &su_m, &lcc_a_m, local_global_index,
				&nsend_a, nsend, sendl, epart_send, epart_recv, &sendl_m, &epart_send_m, &epart_recv_m,
				&l_g, &displs, &displs2, &nintcf_a2, nprocs);

	initVar(nintci, nintcf, nextci, nextcf, lcc, &lcc_a,
			*cnorm, *var, *cgup, *bp, *bs, *be, *bn, *bw, *bh, *bl);

	createList(nsend, *nextcf - *nintcf, &epart_send, &epart_recv, &sendl, nprocs, *nintcf + 1,
			nghb_cnt, nghb_to_rank, send_cnt, recv_cnt, send_lst, recv_lst);

	return 0;
}

