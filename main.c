#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define RANDOM_SEED 1010

#define STARTING_DIMENSION 1024
#define MAX_DIMENSION 1024
#define TEST_EVERY_N_TH 128

#define REPEAT 1

double** getMatrix(int option, int dimension)
{
	double** matrix;

	// allocate memory
	matrix = (double**)malloc(dimension * sizeof(double*));

	for(int i=0; i<dimension; i++)
		matrix[i] = (double*)malloc(dimension * sizeof(double));

	// initialise entries
	for(int i=0; i<dimension; i++)
	{
		for(int j=0; j<dimension; j++)
		{
			switch(option)
			{
				case 0:
					matrix[i][j] = 0;
				case 1:
					matrix[i][j] = (rand() % 100 - 50) / 10;
			}
		}
	}

	return matrix;
}

void freeMatrix(double **matrix, int dimension)
{
	for(int i=0; i<dimension; i++)
		free(matrix[i]);

	free(matrix);
}

void transpose(double **matrix, int dimension)
{
	for(int i=0; i<dimension; i++)
	{
		for(int j=0; j<dimension; j++)
		{
			if(i < j)
			{
				double temp = matrix[i][j];
				matrix[i][j] = matrix[j][i];
				matrix[j][i] = temp;
			}
		}
	}
}

void printMatrix(double **matrix)
{
	for(int i=0; i<MAX_DIMENSION; i++)
	{
		for(int j=0; j<MAX_DIMENSION; j++)
			printf("%5.1f ", matrix[i][j]);

		printf("\n");
	}

	printf("\n");
}

void dot(double **A, double **B, double **C, int i, int j, int k)
{
	register double c0=0, c1=0, c2=0, c3=0;

	register double a0 = A[k][i];
	c0 += a0 * B[k][j];
	c1 += a0 * B[k][j+1];
	c2 += a0 * B[k][j+2];
	c3 += a0 * B[k][j+3];

	register double a1 = A[k+1][i];
	c0 += a1 * B[k+1][j];
	c1 += a1 * B[k+1][j+1];
	c2 += a1 * B[k+1][j+2];
	c3 += a1 * B[k+1][j+3];

	register double a2 = A[k+2][i];
	c0 += a2 * B[k+2][j];
	c1 += a2 * B[k+2][j+1];
	c2 += a2 * B[k+2][j+2];
	c3 += a2 * B[k+2][j+3];

	register double a3 = A[k+3][i];
	c0 += a3 * B[k+3][j];
	c1 += a3 * B[k+3][j+1];
	c2 += a3 * B[k+3][j+2];
	c3 += a3 * B[k+3][j+3];

	C[i][j] += c0;
	C[i][j+1] += c1;
	C[i][j+2] += c2;
	C[i][j+3] += c3;
}

double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);
	transpose(A, dimension);

	for(int i=0; i<dimension; i+=4)
	{
		for(int k=0; k<dimension; k+=4)
		{
			for(int j=0; j<dimension; j+=4)
			{
				dot(A, B, C, i, j, k);
				dot(A, B, C, i+1, j, k);
				dot(A, B, C, i+2, j, k);
				dot(A, B, C, i+3, j, k);
			}
		}
	}

	return C;
}

double** gemmm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int j=0; j<dimension; j++)
		{
			for(int k=0; k<dimension; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return C;
}

int main()
{
	// FILE *file = fopen("results.txt", "wt");
	// fclose(file);

	for(int dimension=STARTING_DIMENSION; dimension<=MAX_DIMENSION; dimension+=TEST_EVERY_N_TH)
	{
		// Initialisation
		srand(RANDOM_SEED);
		double **matrixA = getMatrix(1, dimension);
		double **matrixB = getMatrix(1, dimension);
		double **matrixC;

		// Run GeMM and measure the time
		double timesMeasured[REPEAT];

		for(int take=0; take<REPEAT; take++)
		{
			double tStart = clock();
			matrixC = gemm(matrixA, matrixB, dimension);
			double tEnd = clock();

			double timeTaken = (tEnd - tStart) / CLOCKS_PER_SEC;
			timesMeasured[take] = timeTaken;
		}

		// Calculate Avg and SD
		double avg = 0, sd = 0;

		for(int i=0; i<REPEAT; i++)
			avg += timesMeasured[i];

		avg /= REPEAT;

		for(int i=0; i<REPEAT; i++)
			sd += pow(timesMeasured[i] - avg, 2);

		sd = sqrt(sd / REPEAT);

		printf("size=%3d, %d tries, avg = %.4f, SD = %.4f\n", dimension, REPEAT, avg, sd);

		// Save to file
		FILE *file = fopen("results.txt", "at");
		fprintf(file, "size = %d, avg time = %.4f, SD = %.4f\n", dimension, avg, sd);
		fclose(file);

		// Post-initialisation
		// printMatrix(matrixA);
		// printMatrix(matrixB);
		// printMatrix(matrixC);

		freeMatrix(matrixA, dimension);
		freeMatrix(matrixB, dimension);
		freeMatrix(matrixC, dimension);
	}
	
	return 0;
}