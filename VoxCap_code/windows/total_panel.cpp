
#include <string.h>
#include <iostream>
//#include <algorithm>
//#include <vector>
//#include <malloc.h>
#include "matrix.h"
#include "mex.h"
using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
	double* boolean_tens = (double*)mxGetPr(prhs[0]);
	double* grid_intcon = (double*)mxGetPr(prhs[1]);
	double dx = mxGetScalar(prhs[2]);

	const mwSize *dim_array;
	dim_array = mxGetDimensions(prhs[0]);
	int L, M, N;
	L = *dim_array;
	M = *(dim_array + 1);
	N = *(dim_array + 2);
	mwSize four_D[4] = { L,M,N,24 };
	mwSize *dims;
	dims = four_D;

	plhs[0] = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
	double* total_panels = (double*)mxGetData(plhs[0]);


	int mm, ll, kk;

	for (kk = 0; kk < L; kk++)
	{
		for (ll = 0; ll < M; ll++)
		{
			for (mm = 0; mm < N; mm++)
			{
			
				if (std::abs(boolean_tens[L*M*mm + ll * L + kk] -1.0) < 1e-12)
				{
					total_panels[L*M*N * 0 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 0 + L * M*mm + ll * L + kk] - dx / 2.0;
					total_panels[L*M*N * 1 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 1 + L * M*mm + ll * L + kk];
					total_panels[L*M*N * 2 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 2 + L * M*mm + ll * L + kk];
					total_panels[L*M*N * 3 + L * M*mm + ll * L + kk] = 1;
					total_panels[L*M*N * 8 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 0 + L * M*mm + ll * L + kk];
					total_panels[L*M*N * 9 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 1 + L * M*mm + ll * L + kk] - dx / 2.0;
					total_panels[L*M*N * 10 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 2 + L * M*mm + ll * L + kk];
					total_panels[L*M*N * 11 + L * M*mm + ll * L + kk] = 2;
					total_panels[L*M*N * 16 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 0 + L * M*mm + ll * L + kk];
					total_panels[L*M*N * 17 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 1 + L * M*mm + ll * L + kk];
					total_panels[L*M*N * 18 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 2 + L * M*mm + ll * L + kk] - dx / 2.0;
					total_panels[L*M*N * 19 + L * M*mm + ll * L + kk] = 3;

					if (((kk + 1) <= (L - 1) && boolean_tens[L*M*mm + ll * L + (kk+1)] < 1e-12) || kk == (L - 1))
					{
						total_panels[L*M*N * 4 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 0 + L * M*mm + ll * L + kk] + dx / 2.0;
						total_panels[L*M*N * 5 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 1 + L * M*mm + ll * L + kk];
						total_panels[L*M*N * 6 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 2 + L * M*mm + ll * L + kk];
						total_panels[L*M*N * 7 + L * M*mm + ll * L + kk] = 1;
					
					}
					if (((ll + 1) <= (M - 1) && boolean_tens[L*M*mm + (ll+1) * L + kk] < 1e-12) || ll == (M - 1))
					{
						total_panels[L*M*N * 12 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 0 + L * M*mm + ll * L + kk];
						total_panels[L*M*N * 13 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 1 + L * M*mm + ll * L + kk] + dx / 2.0;
						total_panels[L*M*N * 14 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 2 + L * M*mm + ll * L + kk];
						total_panels[L*M*N * 15 + L * M*mm + ll * L + kk] = 2;
						
					}
					if (((mm + 1) <= (N - 1) && boolean_tens[L*M*(mm+1) + ll * L + kk] < 1e-12) || mm ==(N - 1))
					{
						total_panels[L*M*N * 20 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 0 + L * M*mm + ll * L + kk];
						total_panels[L*M*N * 21 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 1 + L * M*mm + ll * L + kk];
						total_panels[L*M*N * 22 + L * M*mm + ll * L + kk] = grid_intcon[L*M*N * 2 + L * M*mm + ll * L + kk] + dx / 2.0;
						total_panels[L*M*N * 23 + L * M*mm + ll * L + kk] = 3;
					}

				}
			}
		}
	}
}

//L*M*N * 0 + L * M*mm + ll * L + kk
//L*M*mm + kk * M + ll
//L*M*N * 0 + L * M*mm + kk * M + ll


//M*N * 24 * kk + N * 24 * ll + 24 * mm + 18
//M * N * kk + N * ll + mm
//M*N * 3 * kk + N * 3 * ll + 3 * mm + 2






