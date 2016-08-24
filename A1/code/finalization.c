/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include <string.h>
#include "util_write_files.h"
#include "vol2mesh.h"
#include "vol2mesh.c"

void finalization(char* file_in, char* prefix, int total_iters,
		double residual_ratio, int nintci, int nintcf, double* var,
		double* cgup, double* su, int **lcc, float ptime, float mflops) {

	int nodeCnt;
	int **points, **elems;
	char file_out[30] = "summary.out";
	
	int status = write_result(file_in, file_out, nintci, nintcf, var,
			total_iters, residual_ratio, ptime, mflops);

	if (status != 0)
		fprintf(stderr, "Error when trying to write to file %s\n", file_out);

	vol2mesh(nintci, nintcf, lcc, &nodeCnt, &points, &elems);

	char file_su[20];
	char file_var[20];
	char file_cgup[20];
	
	strcpy(file_su, prefix);
	strcat(file_su, ".SU.vtk");
//	printf("%s",file_su);
	strcpy(file_var, prefix);
	strcat(file_var, ".VAR.vtk");
//	printf("%s",file_su);
	strcpy(file_cgup, prefix);
	strcat(file_cgup, ".CGUP.vtk");
//	printf("%s",file_su);
	
	printf("\nWriting vtk file -- %s\n",file_su);
	write_result_vtk(file_su, nintci, nintcf, nodeCnt, points, elems, su);
	printf("Writing vtk file -- %s\n",file_var);
	write_result_vtk(file_var, nintci, nintcf, nodeCnt, points, elems, var);
	printf("Writing vtk file -- %s\n\n",file_cgup);
	write_result_vtk(file_cgup, nintci, nintcf, nodeCnt, points, elems, cgup);
}

