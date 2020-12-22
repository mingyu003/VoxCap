
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

	int L= mxGetScalar(prhs[0]);
	int M = mxGetScalar(prhs[1]);
	int N = mxGetScalar(prhs[2]);

	mwSize three_D_x[3] = { L+1,M,N };
	mwSize *dims_x;
	dims_x = three_D_x;
	mwSize three_D_y[3] = { L,M+1,N };
	mwSize *dims_y;
	dims_y = three_D_y;
	mwSize three_D_z[3] = { L,M,N+1 };
	mwSize *dims_z;
	dims_z = three_D_z;


	plhs[0] = mxCreateNumericArray(3, dims_x, mxDOUBLE_CLASS, mxREAL);
	plhs[1] = mxCreateNumericArray(3, dims_y, mxDOUBLE_CLASS, mxREAL);
	plhs[2] = mxCreateNumericArray(3, dims_z, mxDOUBLE_CLASS, mxREAL);

	double* tensor_ijk_xdir_pnl = (double*)mxGetData(plhs[0]);
	double* tensor_ijk_ydir_pnl = (double*)mxGetData(plhs[1]);
	double* tensor_ijk_zdir_pnl = (double*)mxGetData(plhs[2]);
	int dum_cnt = 1;
	for (int mm = 0; mm<N; mm++)
	{
		for (int ll = 0; ll < M; ll++)
		{
			for (int kk = 0; kk < L+1; kk++)
			{
				tensor_ijk_xdir_pnl[(L+1)*M*mm+(L+1)*ll+kk] = dum_cnt;
				dum_cnt += 1;
			}
		}
	}

	for (int mm = 0; mm < N; mm++)
	{
		for (int ll = 0; ll < M+1; ll++)
		{
			for (int kk = 0; kk < L; kk++)
			{
				tensor_ijk_ydir_pnl[L*(M+1)*mm + L * ll + kk] = dum_cnt;
				dum_cnt += 1;
			}
		}
	}
	for (int mm = 0; mm < N+1; mm++)
	{
		for (int ll = 0; ll < M; ll++)
		{
			for (int kk = 0; kk < L; kk++)
			{
				tensor_ijk_zdir_pnl[L*M*mm + L * ll + kk] = dum_cnt;
				dum_cnt += 1;
			}
		}
	}
	
}














