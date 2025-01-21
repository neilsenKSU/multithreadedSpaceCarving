/*=================================================================
*
* CarveItPosixThreads.c	Simple threaded space carve
*
* Usage: ./CarveItPosixThreads inputFile outputFile [n]
*
*        m = number of tasks = number of masks and projections in input, max 36
*        n = number of threads, default 36
*
* Input file format:
* int sX - voxel size in x direction (100 for seed app)
* int sY - voxel size in y direction (100 for seed app)
* int sZ - voxel size in z direction (100 for seed app)
* unsigned char *V - voxels (each voxel one byte, sX*sY*sZ bytes)
*   -- initially, all 1's, so this could be initialized internally
* int m - number of images (36 as default for seed app)
* repeat m times:
*   double *P - projection matrix (12 doubles each) 
*   int smX - mask size in x direction (200 for seed app)
*   int smY - mask size in y direction (200 for seed app)
*   unsigned char *mask - mask matrix (smX * smY bytes)
*   -- note that masks are all assumed to be the same size for storage/usage
* double width - physical width (10.0 mm for seed app)
* double height - physical height (10.0 mm for seed app)
* double depth - physical depth (10.0 mm for seed app)
*=================================================================*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define MAX_M 36

unsigned int sX, sY, sZ;
unsigned char *voxels;
double *P_ptr[36];
unsigned int smX, smY;
unsigned char *mask_ptr[36];
double width, height, depth;
pthread_t tid[36];
int n = 36; // default to 36 threads
int m; // number of tasks = number of masks and projection in input, max 36
       // ASSUME n <= m
       // thread 0 executes task 0, n, 2*n, .. <= m
       // e.g., with 4 threads, each has 9 tasks if m = 36

// void carve(
//	unsigned int sX, unsigned int sY, unsigned int sZ,
//	unsigned char *V,
//	double *P,
//	unsigned int smX, unsigned int smY,
//	unsigned char *mask,
//	double width, double height, double depth
//	)

void carve(uint64_t i)
{
	char *V = voxels;
	double *P = P_ptr[i];
	unsigned char *mask = mask_ptr[i];
	double EPS = 1e-12;
	double dx = width / ((double)(sX));
	double dy = height / ((double)(sY));
	double dz = depth / ((double)(sZ));
	double x0 = -width / 2 + dx / 2.0;
	double y0 = -height / 2 + dy / 2.0;
	double z0 = -depth / 2 + dz / 2.0;

	int iqx, iqy;

	unsigned int iz = 0, ix = 0, iy = 0;

	for (iz = 0; iz<sZ; iz++){
		for (ix = 0; ix<sX; ix++){
			for (iy = 0; iy<sY; iy++){
				double x = x0 + (double)(ix + 1)*dx;
				double y = y0 + (double)(iy + 1)*dy;
				double z = z0 + (double)(iz + 1)*dz;

				double qx = P[0 + 3 * 0] * x + P[0 + 3 * 1] * y + P[0 + 3 * 2] * z + P[0 + 3 * 3];
				double qy = P[1 + 3 * 0] * x + P[1 + 3 * 1] * y + P[1 + 3 * 2] * z + P[1 + 3 * 3];
				double qz = P[2 + 3 * 0] * x + P[2 + 3 * 1] * y + P[2 + 3 * 2] * z + P[2 + 3 * 3];

				if (fabs(qz)<EPS){
					continue;
				}
				iqx = (int)floor(qx / qz + 0.5);
				iqy = (int)floor(qy / qz + 0.5);

				if (iqx >= 0 && iqx<smX && iqy >= 0 && iqy<smY){
					if (mask[iqy + smY*iqx] == 0){
						V[iy + sY*(ix + sX*iz)] = 0;
					}
				}
				else {
					V[iy + sY*(ix + sX*iz)] = 0;
				}
			}
		}
	}
	return;
}

void *carveThread(void *arg)
{
	uint64_t i = (uint64_t) arg;
	carve(i);
	i+=n;
	while(i<m) {
		carve(i);
		i+=n;
	}
	return NULL;
}

//Input file layout binary without spaces:
//sX sY sZ Voxels P smX smY mask width height depth
void inputData(char* i_filename,
	unsigned int* o_sX, unsigned int* o_sY, unsigned int* o_sZ,
	unsigned char** o_voxels,
	double** o_P,
	unsigned int* o_smX, unsigned int* o_smY,
	unsigned char** o_mask,
	double* o_width, double* o_height, double* o_depth, int* m)
{
	int i, mm, n;
	FILE *fp;
	fp = fopen(i_filename, "rb");

	fread(o_sX, sizeof(unsigned int), 1, fp);
	fread(o_sY, sizeof(unsigned int), 1, fp);
	fread(o_sZ, sizeof(unsigned int), 1, fp);

	*o_voxels = (unsigned char*)malloc(sizeof(unsigned char)* (*o_sX) * (*o_sY) * (*o_sZ));
	fread(*o_voxels, sizeof(unsigned char), (*o_sX) * (*o_sY) * (*o_sZ), fp);
        fread(m, sizeof(int), 1, fp);
        mm = (*m);
	for (i=0; i<mm; i++)
	{
		if (i==0)
			*o_P = (double*)malloc(sizeof(double)*3*4*mm); 
		fread((*o_P+12*i), sizeof(double), 3 * 4, fp);
		fread(o_smX, sizeof(unsigned int), 1, fp);
		fread(o_smY, sizeof(unsigned int), 1, fp);
		n = (*o_smX)*(*o_smY);
		if (i==0)
			*o_mask = (unsigned char*)malloc(sizeof(unsigned char)*n*mm);
		fread((*o_mask+n*i), sizeof(unsigned char), n, fp);
	}
	fread(o_width, sizeof(double), 1, fp);
	fread(o_height, sizeof(double), 1, fp);
	fread(o_depth, sizeof(double), 1, fp);

	fclose(fp);
}

//Output file layout binary without spaces:
//Voxels
void outputData(char* i_filename, unsigned int n_voxels, unsigned char *voxels)
{
	FILE *fp;
	fp = fopen(i_filename, "wb");

	fwrite(voxels, sizeof(unsigned char), n_voxels, fp);

	fclose(fp);
}

int main(int argc, char *argv[])
{
	char *inputFileName, *outputFileName;
        char *fileName;
	double *P;
	unsigned char *mask;
	uint64_t i;
	int iz = 0, ix = 0, iy = 0;
	int size;
 
        if (argc > 1)
		inputFileName = argv[1];
        if (argc > 2)
		outputFileName = argv[2];
        if (argc > 3)
		n = atoi(argv[3]);
	else
		n = 36;
	// ASSUME: n <= m

        fileName = malloc(strlen(outputFileName)+3);
	//Read input from input exchange file
	inputData(inputFileName, &sX, &sY, &sZ, &voxels, &P,
		&smX, &smY, &mask, &width, &height, &depth,&m);
 	for (i=0; i<m; i++) {
		P_ptr[i] = (&P[12*i]);
		mask_ptr[i] = (&mask[40000*i]);
	}
 	for (i=0; i<n; i++) {
                pthread_create(&tid[i], NULL, carveThread, (void *) i);
		// carve(sX, sY, sZ, voxels, P_ptr[i], smX, smY, mask_ptr[i], width, height, depth);
		// P_ptr += 12;
		// mask_ptr += 200*200;
		// Added for debugging
		// printf("sX %u, sY %u, sZ %u, voxels %d, P_ptr %d\n",sX,sY,sZ,voxels,P_ptr);
                // printf("smX %u, smY %u, mask_ptr %d \n", smX, smY, mask_ptr);
            	// printf("width %lf, height %lf, depth %lf\n", width, height, depth);
      		// sprintf(fileName,"%s%02d",outputFileName,i);
		// outputData(fileName, sX * sY * sZ, voxels);
	}
        for (i=0; i<n; i++)
	 	pthread_join(tid[i],NULL);

	size = 0;
	for (iz = 0; iz<sZ; iz++)
		for (ix = 0; ix<sX; ix++)
			for (iy = 0; iy<sY; iy++)
				size += (int) voxels[iy + sY*(ix + sX*iz)]; 

	printf("Size = %d\n", size);
	//Write voxel grid to output exchange file
	outputData(outputFileName, sX * sY * sZ, voxels);

	return 0;
}
