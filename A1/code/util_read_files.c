/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "util_read_files.h"

int read_formatted(char *filename, char *format, int *nintci, int *nintcf,
		int *nextci, int *nextcf, int ***lcc, double **bs, double **be,
		double **bn, double **bw, double **bl, double **bh, double **bp,
		double **su) {

	int i;

	if (!strcmp(format, "text")) {
		/*READING VALUES FROM TEXT FILE*/
		FILE *fp = fopen(filename, "r");
		if (fp == NULL) {
			printf("Error opening file '%s'\nMake sure the file is a text file\n", filename);
			return -1;
		}

		//4 variables in total!!!
		fscanf(fp, "%d", nintci);
		fscanf(fp, "%d", nintcf);
		fscanf(fp, "%d", nextci);
		fscanf(fp, "%d", nextcf);

		//allocating lcc
		if ((*lcc = (int**) malloc(((*nintcf) - (*nintci) + 1) * sizeof(int*)))
				== NULL) {
			printf("malloc failed to allocate first dimension of lcc (nintcf)");
			return -1;
		}
		for (i = (*nintci); i <= (*nintcf); i++) {
			if (((*lcc)[i] = (int *) malloc(6 * sizeof(int))) == NULL) {
				printf("malloc(lcc) failed\n");
				return -1;
			}
		}

		//start reading lcc
		//note that c array index starts from 0 while fortran starts from 1!
		for (i = (*nintci); i <= (*nintcf); i++) {
			fscanf(fp, "%d", &(*lcc)[i][0]);
			fscanf(fp, "%d", &(*lcc)[i][1]);
			fscanf(fp, "%d", &(*lcc)[i][2]);
			fscanf(fp, "%d", &(*lcc)[i][3]);
			fscanf(fp, "%d", &(*lcc)[i][4]);
			fscanf(fp, "%d", &(*lcc)[i][5]);
		}

		// allocate other arrays
		if ((*bs = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*be = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bn = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bw = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bl = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bh = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bp = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*su = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}

		// read the other arrays
		for (i = (*nintci); i <= *nintcf; i++) {
			fscanf(fp, "%lf", &((*bs)[i]));
			fscanf(fp, "%lf", &((*be)[i]));
			fscanf(fp, "%lf", &((*bn)[i]));
			fscanf(fp, "%lf", &((*bw)[i]));
			fscanf(fp, "%lf", &((*bl)[i]));
			fscanf(fp, "%lf", &((*bh)[i]));
			fscanf(fp, "%lf", &((*bp)[i]));
			fscanf(fp, "%lf", &((*su)[i]));
		}
		fclose(fp);
		return 0;
	} else if (!strcmp(format, "bin")) {


		/*READING VALUES FROM BINARY FILE*/
		FILE *fp = fopen(filename, "rb");
		if (fp == NULL) {
			printf("Error opening file %s\nMake sure the file is '.bin'\n", filename);
			return -1;
		}

		fread(nintci, sizeof(int), 1, fp);
		fread(nintcf, sizeof(int), 1, fp);
		fread(nextci, sizeof(int), 1, fp);
		fread(nextcf, sizeof(int), 1, fp);

		//allocating lcc
		if ((*lcc = (int**) malloc(((*nintcf) - (*nintci) + 1) * sizeof(int*)))
				== NULL) {
			printf("malloc failed to allocate first dimension of lcc (nintcf)");
			return -1;
		}

		for (i = (*nintci); i <= (*nintcf); i++) {
			if (((*lcc)[i] = (int *) malloc(6 * sizeof(int))) == NULL) {
				printf("malloc(lcc) failed\n");
				return -1;
			}
		}

		for (i = *nintci; i <= *nintcf; i++) {
			fread((*lcc)[i], sizeof(int), 6, fp);
		}
		// allocate other arrays
		if ((*bs = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*be = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bn = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bw = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bl = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bh = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*bp = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		if ((*su = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL) {
			printf("malloc() failed\n");
			return -1;
		}
		fread(*bs, sizeof(double), (*nextcf + 1), fp);
		fread(*be, sizeof(double), (*nextcf + 1), fp);
		fread(*bn, sizeof(double), (*nextcf + 1), fp);
		fread(*bw, sizeof(double), (*nextcf + 1), fp);
		fread(*bl, sizeof(double), (*nextcf + 1), fp);
		fread(*bh, sizeof(double), (*nextcf + 1), fp);
		fread(*bp, sizeof(double), (*nextcf + 1), fp);
		fread(*su, sizeof(double), (*nextcf + 1), fp);

		fclose(fp);
		return 0;

	} else {
		printf("Please Enter bin or text file !!!\n");
		return -1;
	}

}
