/**************************************************
	
	Collection of all the gemm functions tested

***************************************************/

// initial
// size=1024, 5 tries, avg = 4.8129, SD = 0.1563
double** gemm(double **A, double **B, dimension)
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

// swap j and k loops
// size=1024, 5 tries, avg = 2.2263, SD = 0.0067
double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k++)
		{
			for(int j=0; j<dimension; j++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return C;
}

// groups of 1x4
// size=1024, 5 tries, avg = 2.0083, SD = 0.0099
double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k++)
		{
			for(int j=0; j<dimension; j+=4)
			{
				C[i][j] += A[i][k] * B[k][j];
				C[i][j+1] += A[i][k] * B[k][j+1];
				C[i][j+2] += A[i][k] * B[k][j+2];
				C[i][j+3] += A[i][k] * B[k][j+3];
			}
		}
	}

	return C;
}

// using register, groups of 1x4
// size=1024, 5 tries, avg = 1.6724, SD = 0.0069
double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k++)
		{
			for(int j=0; j<dimension; j+=4)
			{
				register double repeatedA = A[i][k];

				C[i][j] += repeatedA * B[k][j];
				C[i][j+1] += repeatedA * B[k][j+1];
				C[i][j+2] += repeatedA * B[k][j+2];
				C[i][j+3] += repeatedA * B[k][j+3];
			}
		}
	}

	return C;
}

// using register, groups of 1x8
// size=1024, 5 tries, avg = 1.5148, SD = 0.0062
double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k++)
		{
			for(int j=0; j<dimension; j+=8)
			{
				register double repeatedA = A[i][k];

				C[i][j] += repeatedA * B[k][j];
				C[i][j+1] += repeatedA * B[k][j+1];
				C[i][j+2] += repeatedA * B[k][j+2];
				C[i][j+3] += repeatedA * B[k][j+3];
				C[i][j+4] += repeatedA * B[k][j+4];
				C[i][j+5] += repeatedA * B[k][j+5];
				C[i][j+6] += repeatedA * B[k][j+6];
				C[i][j+7] += repeatedA * B[k][j+7];
			}
		}
	}

	return C;
}

// with the 'dot' function, groups of 1x4
// size=1024, 1 tries, avg = 1.9722, SD = 0.0000
void dot(double **A, double **B, double **C, int i, int j, int k)
{
	register double repeatedA = A[i][k];

	C[i][j] += repeatedA * B[k][j];
	C[i][j+1] += repeatedA * B[k][j+1];
	C[i][j+2] += repeatedA * B[k][j+2];
	C[i][j+3] += repeatedA * B[k][j+3];
}

double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k++)
		{
			for(int j=0; j<dimension; j+=4)
			{
				dot(A, B, C, i, j, k);
			}
		}
	}

	return C;
}

// 1x4 with 'dot' function
// size=1024, 1 tries, avg = 1.8326, SD = 0.0000
void dot(double **A, double **B, double **C, int i, int j, int kprime)
{
	for(int k=kprime; k<kprime+4; k++)
	{
		register double repeatedA = A[i][k];

		C[i][j] += repeatedA * B[k][j];
		C[i][j+1] += repeatedA * B[k][j+1];
		C[i][j+2] += repeatedA * B[k][j+2];
		C[i][j+3] += repeatedA * B[k][j+3];	
	}
}

double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k+=4)
		{
			for(int j=0; j<dimension; j+=4)
			{
				dot(A, B, C, i, j, k);
			}
		}
	}

	return C;
}

// 1x4 with 'dot', register keyword
// size=1024, 1 tries, avg = 1.6519, SD = 0.0000
void dot(double **A, double **B, double **C, int i, int j, int kprime)
{
	register double c0=0, c1=0, c2=0, c3=0;

	for(int k=kprime; k<kprime+4; k++)
	{
		register double repeatedA = A[i][k];

		c0 += repeatedA * B[k][j];
		c1 += repeatedA * B[k][j+1];
		c2 += repeatedA * B[k][j+2];
		c3 += repeatedA * B[k][j+3];
	}

	C[i][j] += c0;
	C[i][j+1] += c1;
	C[i][j+2] += c2;
	C[i][j+3] += c3;
}

double** gemm(double **A, double **B, int dimension)
{
	double **C = getMatrix(0, dimension);

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k+=4)
		{
			for(int j=0; j<dimension; j+=4)
			{
				dot(A, B, C, i, j, k);
			}
		}
	}

	return C;
}

// 1x4 with 'dot', register keyword, explicit loops
// size=1024, 1 tries, avg = 1.4804, SD = 0.0000
void dot(double **A, double **B, double **C, int i, int j, int k)
{
	register double c0=0, c1=0, c2=0, c3=0;

	register double a0 = A[i][k];
	c0 += a0 * B[k][j];
	c1 += a0 * B[k][j+1];
	c2 += a0 * B[k][j+2];
	c3 += a0 * B[k][j+3];

	register double a1 = A[i][k+1];
	c0 += a1 * B[k+1][j];
	c1 += a1 * B[k+1][j+1];
	c2 += a1 * B[k+1][j+2];
	c3 += a1 * B[k+1][j+3];

	register double a2 = A[i][k+2];
	c0 += a2 * B[k+2][j];
	c1 += a2 * B[k+2][j+1];
	c2 += a2 * B[k+2][j+2];
	c3 += a2 * B[k+2][j+3];

	register double a3 = A[i][k+3];
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

	for(int i=0; i<dimension; i++)
	{
		for(int k=0; k<dimension; k+=4)
		{
			for(int j=0; j<dimension; j+=4)
			{
				dot(A, B, C, i, j, k);
			}
		}
	}

	return C;
}

// 4x4, dot, registers, explicit loops
// size=1024, 1 tries, avg = 1.4364, SD = 0.0000
void dot(double **A, double **B, double **C, int i, int j, int k)
{
	register double c0=0, c1=0, c2=0, c3=0;

	register double a0 = A[i][k];
	c0 += a0 * B[k][j];
	c1 += a0 * B[k][j+1];
	c2 += a0 * B[k][j+2];
	c3 += a0 * B[k][j+3];

	register double a1 = A[i][k+1];
	c0 += a1 * B[k+1][j];
	c1 += a1 * B[k+1][j+1];
	c2 += a1 * B[k+1][j+2];
	c3 += a1 * B[k+1][j+3];

	register double a2 = A[i][k+2];
	c0 += a2 * B[k+2][j];
	c1 += a2 * B[k+2][j+1];
	c2 += a2 * B[k+2][j+2];
	c3 += a2 * B[k+2][j+3];

	register double a3 = A[i][k+3];
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