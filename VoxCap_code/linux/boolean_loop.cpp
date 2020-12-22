
#include <string.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <malloc.h>
#include "matrix.h"
#include "mex.h"
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{

	double* x_bnd = (double*)mxGetPr(prhs[0]);
	double* y_bnd = (double*)mxGetPr(prhs[1]);
	double* z_bnd = (double*)mxGetPr(prhs[2]);
	double* r = (double*)mxGetPr(prhs[3]);

	const mwSize *dim_array;

	dim_array = mxGetDimensions(prhs[3]);
	mwSize L, M, N;
	L = *dim_array;
	M = *(dim_array + 1);
	N = *(dim_array + 2);

	mwSize three_D[3] = { L,M,N };
	mwSize *dims;
	dims = three_D;
	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);

	double* boolean_tens = (double*)mxGetData(plhs[0]);

	double tola = 1e-12;

	for (int ll = 0; ll < L; ll++)
	{
		for (int mm = 0; mm < M; mm++)
		{
			for (int nn = 0; nn < N; nn++)
			{
				if (r[(L*M*N) * 0 + (L*M)*nn + mm * L + ll] > x_bnd[0] - tola && r[(L*M*N) * 0 + (L*M)*nn + mm * L + ll]< x_bnd[1] + tola && \
					r[(L*M*N) * 1 + (L*M)*nn + mm * L + ll]> y_bnd[0] - tola && r[(L*M*N) * 1 + (L*M)*nn + mm * L + ll] < y_bnd[1] + tola && \
					r[(L*M*N) * 2 + (L*M)*nn + mm * L + ll] > z_bnd[0] - tola && r[(L*M*N) * 2 + (L*M)*nn + mm * L + ll] < z_bnd[1] + tola)
				{
					boolean_tens[(L*M)*nn + mm * L + ll] = 1;
				}
			}
		}
	}
}


