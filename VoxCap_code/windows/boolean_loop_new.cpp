
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
    double dx = mxGetScalar(prhs[4]);

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

    double org_comp_domain[3];
    
    mwSize x_low, x_up, y_low, y_up, z_low, z_up;
    if (round((x_bnd[0]-tola-org_comp_domain[0])/dx)< 0)
    {
        x_low=0;
    }
    else
    {
        x_low=round((x_bnd[0]-tola-org_comp_domain[0])/dx);
    }
    if (round((x_bnd[1]+tola+org_comp_domain[0])/dx) > L)
    {
        x_up=L;
    }
    else
    {
        x_up=round((x_bnd[1]+tola+org_comp_domain[0])/dx);
    }       
    if (round((y_bnd[0]-tola-org_comp_domain[1])/dx)<0)
    {
        y_low=0;
    }
    else
    {
        y_low=round((y_bnd[0]-tola-org_comp_domain[1])/dx);
    }            
    if (round((y_bnd[1]+tola+org_comp_domain[1])/dx) > M)
    {
        y_up=M;
    }
    else
    {
        y_up=round((y_bnd[1]+tola+org_comp_domain[1])/dx);
    } 
    if (round((z_bnd[0]-tola-org_comp_domain[2])/dx)< 0)
    {
        z_low=0;
    }
    else
    {
        z_low=round((z_bnd[0]-tola-org_comp_domain[2])/dx);
    } 
    if (round((z_bnd[1]+tola+org_comp_domain[2])/dx) > N)
    {
        z_up=N;
    }
    else
    {
        z_up=round((z_bnd[1]+tola+org_comp_domain[2])/dx);
    } 

    org_comp_domain[0]=r[(L*M*N) * 0 + (L*M)*0 + 0 * L + 0];
    org_comp_domain[1]=r[(L*M*N) * 1 + (L*M)*0 + 0 * L + 0];
    org_comp_domain[2]=r[(L*M*N) * 2 + (L*M)*0 + 0 * L + 0];
   
	for (int ll = x_low; ll < x_up; ll++)
	{
		for (int mm = y_low; mm < y_up; mm++)
		{
			for (int nn = z_low; nn < z_up; nn++)
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


