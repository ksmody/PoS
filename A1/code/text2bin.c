#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{

int i;
char *Text, *Bin;
int nintci, nintcf, nextci, nextcf;
int **lcc;
double *bs, *be, *bn, *bw, *bl, *bh, *bp, *su;

//Check for total number of arguments to be 3
if (argc != 3) {
		printf("Please Enter -- binconv <Input file name.dat> <Output file name>\n");
		return -1;
	}

Text = argv[1];
Bin = argv[2];


//Code from UTIL_READ_FILES.c ----------------------
FILE *fp = fopen(Text, "r");
if (fp == NULL)
{
printf("Error opening file %s\n", Text);
return -1;
}
//4 variables in total!!!
fscanf(fp, "%d", &nintci);
fscanf(fp, "%d", &nintcf);
fscanf(fp, "%d", &nextci);
fscanf(fp, "%d", &nextcf);
//allocating lcc
if ((lcc = (int**) malloc( ( (nintcf) - (nintci) + 1) * sizeof(int*) )) == NULL)
{
printf("malloc failed to allocate first dimension of lcc (nintcf)");
return -1;
}
for (i = (nintci); i <= (nintcf); i++)
{
if ((lcc[i] = (int *) malloc( 6 * sizeof(int) )) == NULL)
{
printf("malloc(lcc) failed\n");
return -1;
}
}
//start reading lcc
//note that c array index starts from 0 while fortran starts from 1!
for (i = nintci; i <= nintcf; i++)
{
fscanf(fp, "%d", &(lcc[i][0]));
fscanf(fp, "%d", &(lcc[i][1]));
fscanf(fp, "%d", &(lcc[i][2]));
fscanf(fp, "%d", &(lcc[i][3]));
fscanf(fp, "%d", &(lcc[i][4]));
fscanf(fp, "%d", &(lcc[i][5]));
}
// allocate other arrays
if ((bs = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((be = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((bn = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((bw = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((bl = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((bh = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((bp = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
if ((su = (double *) malloc((nextcf + 1) * sizeof(double))) == NULL)
{
printf("malloc() failed\n");
return -1;
}
// read the other arrays
for (i = nintci; i <= nintcf; i++)
{
fscanf(fp, "%lf", &(bs[i]));
fscanf(fp, "%lf", &(be[i]));
fscanf(fp, "%lf", &(bn[i]));
fscanf(fp, "%lf", &(bw[i]));
fscanf(fp, "%lf", &(bl[i]));
fscanf(fp, "%lf", &(bh[i]));
fscanf(fp, "%lf", &(bp[i]));
fscanf(fp, "%lf", &(su[i]));
}
fclose(fp);


// Text to Binary Code ------------------------
FILE *fBin = fopen(strcat(Bin, ".bin"), "wb+");
if (!fBin)
{
	printf("Error Opening File %s\n", Bin);
	return -1;
}

fwrite(&nintci, sizeof(int), 1, fBin);
fwrite(&nintcf, sizeof(int), 1, fBin);
fwrite(&nextci, sizeof(int), 1, fBin);
fwrite(&nextcf, sizeof(int), 1, fBin);

for (i = nintci; i <= nintcf; i++)
	{
	fwrite(lcc[i], sizeof(int), 6, fBin);
	}

fwrite(bs, sizeof(double), (nextcf + 1), fBin);
fwrite(be, sizeof(double), (nextcf + 1), fBin);
fwrite(bn, sizeof(double), (nextcf + 1), fBin);
fwrite(bw, sizeof(double), (nextcf + 1), fBin);
fwrite(bl, sizeof(double), (nextcf + 1), fBin);
fwrite(bh, sizeof(double), (nextcf + 1), fBin);
fwrite(bp, sizeof(double), (nextcf + 1), fBin);
fwrite(su, sizeof(double), (nextcf + 1), fBin);

fclose(fBin);


//Reading Binary file to verify Content--------------------------
fBin = fopen(Bin, "rb");
if (!fBin)
{
	printf("Error Opening File %s\n", Bin);
	return -1;
}

fread(&nintci, sizeof(int), 1, fBin);
fread(&nintcf, sizeof(int), 1, fBin);
fread(&nextci, sizeof(int), 1, fBin);
fread(&nextcf, sizeof(int), 1, fBin);

for (i = nintci; i <= nintcf; i++)
	{
	fread(lcc[i], sizeof(int), 6, fBin);
	}

fread(bs, sizeof(double), (nextcf + 1), fBin);
fread(be, sizeof(double), (nextcf + 1), fBin);
fread(bn, sizeof(double), (nextcf + 1), fBin);
fread(bw, sizeof(double), (nextcf + 1), fBin);
fread(bl, sizeof(double), (nextcf + 1), fBin);
fread(bh, sizeof(double), (nextcf + 1), fBin);
fread(bp, sizeof(double), (nextcf + 1), fBin);
fread(su, sizeof(double), (nextcf + 1), fBin);

fclose(fBin);

for (i=0; i <= nintcf - nintci; i++)
	{
	free(lcc[i]);
	}
free(lcc);

free(su);
free(bp);
free(bh);
free(bl);
free(bw);
free(bn);
free(be);
free(bs);
}
