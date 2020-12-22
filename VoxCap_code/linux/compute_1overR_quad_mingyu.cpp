#include <string>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <malloc.h>
#include "compute_1overR_quad_mingyu.h"
#include "mex.h"
#include <iterator>

using namespace std;
template<class T>
int length(T& arr)
{
	return sizeof(arr) / sizeof(arr[0]);
}


template<class T>
double norm(T& arr)
{
	return sqrt(pow(arr[0], 2) + pow(arr[1], 2) + pow(arr[2], 2));
}

double compute_1overR_quad_mingyu(double src_cen[3], double obs_cen[3], double unit_normal_src[3], double unit_normal_obs[3], double smpl_wghts[], double smpl_pnts_temp[], int len_smpl)
{
	double y = 0.0;
	double int_res = 0.0;
	int plane_src = 0; int plane_obs = 0; int form_type = 0; int ii;
	//int len_smpl = length(smpl_wghts);

	double arr_smpl[4096 * 4];
	double* smpl_pnts = &arr_smpl[0];
	for (int ii = 0; ii < len_smpl * 4; ii++)
	{
		arr_smpl[ii] = *(smpl_pnts_temp + ii);
	}


	if (norm(unit_normal_src) - 1.0 > 1e-13)
	{
		y = norm(unit_normal_src);
		for (plane_src = 0; plane_src < 3; plane_src++)
		{
			unit_normal_src[plane_src] /= y;
		}
	}
	if (norm(unit_normal_obs) - 1.0 > 1e-13)
	{
		y = norm(unit_normal_obs);
		for (plane_obs = 0; plane_obs < 3; plane_obs++)
		{
			unit_normal_obs[plane_obs] /= y;
		}
	}

	if ((std::abs(std::abs(unit_normal_src[2]) - 1.0) < 1.0E-13) && (std::abs
	(unit_normal_src[1]) < 1.0E-13) && (std::abs(unit_normal_src[0]) < 1.0E-13))
	{
		plane_src = 1;
	}
	else if ((std::abs(std::abs(unit_normal_src[1]) - 1.0) < 1.0E-13) && (std::
		abs(unit_normal_src[2]) < 1.0E-13) && (std::abs(unit_normal_src[0]) < 1.0E-13))
	{
		plane_src = 2;
	}
	else if ((std::abs(std::abs(unit_normal_src[0]) - 1.0) < 1.0E-13) && (std::abs
	(unit_normal_src[1]) < 1.0E-13) && (std::abs(unit_normal_src[2]) < 1.0E-13))
	{
		plane_src = 3;
	}

	if ((std::abs(std::abs(unit_normal_obs[2]) - 1.0) < 1.0E-13) && (std::abs
	(unit_normal_obs[1]) < 1.0E-13) && (std::abs(unit_normal_obs[0]) < 1.0E-13))
	{
		plane_obs = 1;
	}
	else if ((std::abs(std::abs(unit_normal_obs[1]) - 1.0) < 1.0E-13) && (std::
		abs(unit_normal_obs[2]) < 1.0E-13) && (std::abs(unit_normal_obs[0]) < 1.0E-13))
	{
		plane_obs = 2;
	}
	else if ((std::abs(std::abs(unit_normal_obs[0]) - 1.0) < 1.0E-13) && (std::abs
	(unit_normal_obs[1]) < 1.0E-13) && (std::abs(unit_normal_obs[2]) < 1.0E-13))
	{
		plane_obs = 3;
	}

	if (plane_obs == plane_src)
	{
		form_type = 1;
	}
	else
	{
		form_type = 2;
	}

	if (form_type == 1)
	{
		if (plane_src == 1)
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[0];
				smpl_pnts[len_smpl + ii] += src_cen[0];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[1];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[1];

				y = 1.0 / sqrt(pow(smpl_pnts[ii] - smpl_pnts[len_smpl + ii], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - smpl_pnts[len_smpl * 3 + ii], 2) + pow(obs_cen[2] - src_cen[2], 2));
				y *= smpl_wghts[ii];
				int_res += y;

			}
		}
		else if (plane_src == 2)
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[0];
				smpl_pnts[len_smpl + ii] += src_cen[0];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[2];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[2];
				y = 1.0 / sqrt(pow(smpl_pnts[ii] - smpl_pnts[len_smpl + ii], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - smpl_pnts[len_smpl * 3 + ii], 2) + pow(obs_cen[1] - src_cen[1], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
		else if (plane_src == 3)
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[1];
				smpl_pnts[len_smpl + ii] += src_cen[1];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[2];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[2];
				y = 1.0 / sqrt(pow(smpl_pnts[ii] - smpl_pnts[len_smpl + ii], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - smpl_pnts[len_smpl * 3 + ii], 2) + pow(obs_cen[0] - src_cen[0], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
	}
	else if (form_type == 2)
	{
		if ((plane_src == 1) && (plane_obs == 2))
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[0];
				smpl_pnts[len_smpl + ii] += src_cen[0];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[1];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[2];
				y = 1.0 / sqrt(pow(smpl_pnts[ii] - smpl_pnts[len_smpl + ii], 2) + pow(obs_cen[1] - smpl_pnts[len_smpl * 3 + ii], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - src_cen[2], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
		else if ((plane_src == 1) && (plane_obs == 3))
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += src_cen[0];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[1];
				smpl_pnts[len_smpl + ii] += src_cen[1];
				smpl_pnts[len_smpl * 3 + ii] += obs_cen[2];
				y = 1.0 / sqrt(pow(obs_cen[0] - smpl_pnts[ii], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - smpl_pnts[len_smpl + ii], 2) + pow(smpl_pnts[len_smpl * 3 + ii] - src_cen[2], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
		else if ((plane_src == 2) && (plane_obs == 3))
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += src_cen[0];
				smpl_pnts[len_smpl + ii] += obs_cen[1];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[2];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[2];
				y = 1.0 / sqrt(pow(obs_cen[0] - smpl_pnts[ii], 2) + pow(smpl_pnts[len_smpl + ii] - src_cen[1], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - smpl_pnts[len_smpl * 3 + ii], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
		else if ((plane_src == 2) && (plane_obs == 1))
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[0];
				smpl_pnts[len_smpl + ii] += src_cen[0];
				smpl_pnts[len_smpl * 2 + ii] += obs_cen[1];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[2];
				y = 1.0 / sqrt(pow(smpl_pnts[ii] - smpl_pnts[len_smpl + ii], 2) + pow(smpl_pnts[len_smpl * 2 + ii] - src_cen[1], 2) + pow(obs_cen[2] - smpl_pnts[len_smpl * 3 + ii], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
		else if ((plane_src == 3) && (plane_obs == 1))
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[0];
				smpl_pnts[len_smpl + ii] += obs_cen[1];
				smpl_pnts[len_smpl * 2 + ii] += src_cen[1];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[2];
				y = 1.0 / sqrt(pow(smpl_pnts[ii] - src_cen[0], 2) + pow(smpl_pnts[len_smpl + ii] - smpl_pnts[len_smpl * 2 + ii], 2) + pow(obs_cen[2] - smpl_pnts[len_smpl * 3 + ii], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
		else if ((plane_src == 3) && (plane_obs == 2))
		{
			for (ii = 0; ii < len_smpl; ii++)
			{
				smpl_pnts[ii] += obs_cen[0];
				smpl_pnts[len_smpl * 2 + ii] += src_cen[1];
				smpl_pnts[len_smpl + ii] += obs_cen[2];
				smpl_pnts[len_smpl * 3 + ii] += src_cen[2];
				y = 1.0 / sqrt(pow(smpl_pnts[ii] - src_cen[0], 2) + pow(obs_cen[1] - smpl_pnts[len_smpl * 2 + ii], 2) + pow(smpl_pnts[len_smpl + ii] - smpl_pnts[len_smpl * 3 + ii], 2));
				y *= smpl_wghts[ii];
				int_res += y;
			}
		}
	}
	return int_res;
}

