
#include <string.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <malloc.h>
//#include "compute_1overR_quad_mingyu.h"
//#include "analytical_derivation_1overR.h"
//#include "finite_difference_coefficients.h"
#include "mex.h"
const double eps2 = 1e-37;
const double pi = 4.0*atan(1.0);
template<class T>
double norm(T& arr)
{
	return sqrt(pow(arr[0], 2) + pow(arr[1], 2) + pow(arr[2], 2));
}
using namespace std;
template<class T>
int length(T& arr)
{
	return sizeof(arr) / sizeof(arr[0]);
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

double finite_difference_coefficients(double dx, int num_diff_pnts, double src_cen[3], double obs_cen[3],
	double unit_normal_src[3], double unit_normal_obs[3], double smpl_wghts[], double smpl_pnts[], int len_smpl)
{
	double res_diff;
	double h = 0.1*dx;
	double obs_cen_term1[3]; double obs_cen_term2[3]; double obs_cen_term3[3]; double obs_cen_term4[3]; double obs_cen_term5[3]; double obs_cen_term6[3];
	double term1; double term2; double term3; double term4; double term5; double term6;
	if (num_diff_pnts == 2) //-1,1
	{
		for (int ii = 0; ii < 3; ii++)
		{
			obs_cen_term1[ii] = obs_cen[ii] + unit_normal_obs[ii] * h;
			obs_cen_term2[ii] = obs_cen[ii] - unit_normal_obs[ii] * h;
		}
		term1 = compute_1overR_quad_mingyu(src_cen, obs_cen_term1, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term2 = compute_1overR_quad_mingyu(src_cen, obs_cen_term2, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		res_diff = (term1 - term2) / (2.0 * h);
	}
	else if (num_diff_pnts == 4) //-1,+1,-2,+2
	{
		for (int ii = 0; ii < 3; ii++)
		{
			obs_cen_term1[ii] = obs_cen[ii] - unit_normal_obs[ii] * 2.0* h;
			obs_cen_term2[ii] = obs_cen[ii] - unit_normal_obs[ii] * h;
			obs_cen_term3[ii] = obs_cen[ii] + unit_normal_obs[ii] * h;
			obs_cen_term4[ii] = obs_cen[ii] + unit_normal_obs[ii] * 2.0*h;
		}
		term1 = compute_1overR_quad_mingyu(src_cen, obs_cen_term1, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term2 = compute_1overR_quad_mingyu(src_cen, obs_cen_term2, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term3 = compute_1overR_quad_mingyu(src_cen, obs_cen_term3, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term4 = compute_1overR_quad_mingyu(src_cen, obs_cen_term4, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		res_diff = (term1 - 8.0*term2 + 8.0*term3 - term4) / (12.0 * h);
	}
	else if (num_diff_pnts == 6) //-1,+1,-2,+2,-3,+3
	{
		for (int ii = 0; ii < 3; ii++)
		{
			obs_cen_term1[ii] = obs_cen[ii] - unit_normal_obs[ii] * 3.0 * h;
			obs_cen_term2[ii] = obs_cen[ii] - unit_normal_obs[ii] * 2.0*h;
			obs_cen_term3[ii] = obs_cen[ii] - unit_normal_obs[ii] * h;
			obs_cen_term4[ii] = obs_cen[ii] + unit_normal_obs[ii] * h;
			obs_cen_term5[ii] = obs_cen[ii] + unit_normal_obs[ii] * 2.0*h;
			obs_cen_term6[ii] = obs_cen[ii] + unit_normal_obs[ii] * 3.0 * h;
		}
		term1 = compute_1overR_quad_mingyu(src_cen, obs_cen_term1, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term2 = compute_1overR_quad_mingyu(src_cen, obs_cen_term2, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term3 = compute_1overR_quad_mingyu(src_cen, obs_cen_term3, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term4 = compute_1overR_quad_mingyu(src_cen, obs_cen_term4, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term5 = compute_1overR_quad_mingyu(src_cen, obs_cen_term5, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		term6 = compute_1overR_quad_mingyu(src_cen, obs_cen_term6, unit_normal_src, unit_normal_obs, smpl_wghts, smpl_pnts, len_smpl);
		res_diff = (-term1 + 9.0 * term2 - 45.0 * term3 + 45.0*term4 - 9.0*term5 + term6) / (60.0 * h);
	}

	return res_diff;
}
int sign(double x)
{
	if (x > eps2) { return 1; }
	else if (x < eps2) { return -1; }
	else { return 0; }
}
double analytical_derivation_1overR(double dx, double src_cen[3], double obs_cen[3], double unit_normal_src[3], double unit_normal_obs[3])
{
	

	const double pi = 4.0*atan(1.0);
	double CC_w_sgn; int fl_self_term; int kk; int mm; int ll;
	double term1; double term2; double term3; double term4; double term5; double term6; double term7; double y;
	double term1_1; double term1_2; double term1_3; double term2_1; double term2_2; double term3_1; double term3_2; double term4_1; double term4_2;
	double term7_a; double term7_b; double term7_b1; double term7_b2; double term7_c;
	double int_res = 0.0;
	int fl_norm; int plane_src; int plane_obs; int form_type;
	double a_arr[4]; double bp_arr[4]; double bo_arr[2]; double co_arr[2]; double aij; double bij; double cij; double z;


	if ((std::abs(obs_cen[0] - src_cen[0]) < 1e-13) && (std::abs(obs_cen[1] - src_cen[1]) < 1e-13) && (std::abs(obs_cen[2] - src_cen[2]) < 1e-13)
		&& (std::abs(unit_normal_src[0] - unit_normal_obs[0]) < 1e-13) && (std::abs(unit_normal_src[1] - unit_normal_obs[1]) < 1e-13) && (std::abs(unit_normal_src[2] - unit_normal_obs[2]) < 1e-13))
	{
		fl_self_term = 1;
	}
	else
	{
		fl_self_term = 0;
	}
	if (fl_self_term == 1)
	{
		int_res = 0;//pow(dx, 2) * 2 * pi*((epsa + epsb) / (epsa - epsb));
	}
	if (fl_self_term == 0)
	{
		if (norm(unit_normal_src) - 1.0 > 1e-13)
		{
			y = norm(unit_normal_src);
			for (fl_norm = 0; fl_norm < 3; fl_norm++)
			{
				unit_normal_src[fl_norm] /= y;
			}
		}
		if (norm(unit_normal_obs) - 1.0 > 1e-13)
		{
			y = norm(unit_normal_obs);
			for (fl_norm = 0; fl_norm < 3; fl_norm++)
			{
				unit_normal_obs[fl_norm] /= y;
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
				z = std::abs(obs_cen[2] - src_cen[2]);
				aij = std::abs(obs_cen[0] - src_cen[0]);
				bij = std::abs(obs_cen[1] - src_cen[1]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bp_arr[0] = bij - dx; bp_arr[1] = bij; bp_arr[2] = bij + dx; bp_arr[3] = bij;
				CC_w_sgn = obs_cen[2] - src_cen[2];
			}

			else if (plane_src == 2)
			{
				z = std::abs(obs_cen[1] - src_cen[1]);
				aij = std::abs(obs_cen[2] - src_cen[2]);
				bij = std::abs(obs_cen[0] - src_cen[0]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bp_arr[0] = bij - dx; bp_arr[1] = bij; bp_arr[2] = bij + dx; bp_arr[3] = bij;
				CC_w_sgn = obs_cen[1] - src_cen[1];
			}
			else if (plane_src == 3)
			{
				z = std::abs(obs_cen[0] - src_cen[0]);
				aij = std::abs(obs_cen[1] - src_cen[1]);
				bij = std::abs(obs_cen[2] - src_cen[2]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bp_arr[0] = bij - dx; bp_arr[1] = bij; bp_arr[2] = bij + dx; bp_arr[3] = bij;
				CC_w_sgn = obs_cen[0] - src_cen[0];
			}
			double rkm; double mult_fact;
			for (kk = 0; kk < 4; kk++)
			{
				for (mm = 0; mm < 4; mm++)
				{
					rkm = sqrt(pow(a_arr[kk], 2.0) + pow(bp_arr[mm], 2.0) + pow(z, 2.0));
					mult_fact = pow(-1, (mm + kk));
					//term1
					if ((std::abs(a_arr[kk] + rkm) < eps2) || (std::abs(rkm) < eps2))
					{
						term1_1 = 0;
					}
					else
					{
						term1_1 = 0.5*(pow(bp_arr[mm], 2.0) - pow(z, 2.0))*a_arr[kk] * z / (rkm*(a_arr[kk] + rkm));
					}
					if (abs(a_arr[kk] + rkm) < eps2)
					{
						term1_2 = 0;
					}
					else
					{
						term1_2 = a_arr[kk] * z*std::log(a_arr[kk] + rkm);
					}
					term1 = term1_1 - term1_2;
					//term2
					if ((std::abs(bp_arr[mm] + rkm) < eps2) || std::abs(rkm) < eps2)
					{
						term2_1 = 0;
					}
					else
					{
						term2_1 = 0.5*(pow(a_arr[kk], 2.0) - pow(z, 2.0))*bp_arr[mm] * z / (rkm*(bp_arr[mm] + rkm));
					}
					if (std::abs(bp_arr[mm] + rkm) < eps2)
					{
						term2_2 = 0;
					}
					else
					{
						term2_2 = bp_arr[mm] * z*std::log(bp_arr[mm] + rkm);
					}
					term2 = term2_1 - term2_2;
					//term3
					if (std::abs(rkm) < eps2)
					{
						term3_1 = 0;
					}
					else
					{
						term3_1 = (pow(bp_arr[mm], 2.0) - 2.0 * pow(z, 2.0) + pow(a_arr[kk], 2.0))*z / (6.0 * rkm);
					}
					term3_2 = (2.0 / 3.0)*z*rkm;
					term3 = term3_1 - term3_2;
					//term4
					if (std::abs(rkm*z) < eps2)
					{
						term4_1 = 0;
					}
					else
					{
						term4_1 = a_arr[kk] * bp_arr[mm] * std::atan(a_arr[kk] * bp_arr[mm] / (z*rkm));
					}
					if ((std::abs(rkm) < eps2) || std::abs((pow(z, 2.0)*pow(rkm, 2.0)) + (pow(a_arr[kk], 2.0)*pow(bp_arr[mm], 2.0)) < eps2))
					{
						term4_2 = 0;
					}
					else
					{
						term4_2 = (pow(a_arr[kk], 2.0)*pow(bp_arr[mm], 2.0)*z / rkm)*(pow(z, 2.0) + pow(rkm, 2.0)) / (pow(z, 2.0)*pow(rkm, 2.0) + pow(a_arr[kk], 2.0)*pow(bp_arr[mm], 2.0));
					}
					term4 = term4_1 - term4_2;

					int_res += (mult_fact*(term1 + term2 - term3 - term4));

				}
			}

			if (std::abs((unit_normal_obs[0] + unit_normal_obs[1] + unit_normal_obs[2]) - sign(CC_w_sgn)) > eps2)
			{
				int_res = 0 - int_res;
			}
		}

		else if (form_type == 2)
		{
			if (plane_src == 1 && plane_obs == 2)
			{
				aij = std::abs(obs_cen[0] - src_cen[0]);
				bij = std::abs(obs_cen[2] - src_cen[2]);
				cij = std::abs(obs_cen[1] - src_cen[1]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
				CC_w_sgn = obs_cen[1] - src_cen[1];
			}
			else if (plane_src == 2 && plane_obs == 1)
			{
				aij = std::abs(obs_cen[0] - src_cen[0]);
				bij = std::abs(obs_cen[1] - src_cen[1]);
				cij = std::abs(obs_cen[2] - src_cen[2]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
				CC_w_sgn = obs_cen[2] - src_cen[2];
			}
			else if (plane_src == 3 && plane_obs == 1)
			{
				aij = std::abs(obs_cen[1] - src_cen[1]);
				bij = std::abs(obs_cen[0] - src_cen[0]);
				cij = std::abs(obs_cen[2] - src_cen[2]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
				CC_w_sgn = obs_cen[2] - src_cen[2];
			}
			else if (plane_src == 1 && plane_obs == 3)
			{
				aij = std::abs(obs_cen[1] - src_cen[1]);
				bij = std::abs(obs_cen[2] - src_cen[2]);
				cij = std::abs(obs_cen[0] - src_cen[0]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
				CC_w_sgn = obs_cen[0] - src_cen[0];
			}
			else if (plane_src == 2 && plane_obs == 3)
			{
				aij = std::abs(obs_cen[2] - src_cen[2]);
				bij = std::abs(obs_cen[1] - src_cen[1]);
				cij = std::abs(obs_cen[0] - src_cen[0]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
				CC_w_sgn = obs_cen[0] - src_cen[0];
			}
			else if (plane_src == 3 && plane_obs == 2)
			{
				aij = std::abs(obs_cen[2] - src_cen[2]);
				bij = std::abs(obs_cen[0] - src_cen[0]);
				cij = std::abs(obs_cen[1] - src_cen[1]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
				CC_w_sgn = obs_cen[1] - src_cen[1];
			}
			double rkml; double mult_fact;
			for (kk = 0; kk < 4; kk++)
			{
				for (mm = 0; mm < 2; mm++)
				{
					for (ll = 0; ll < 2; ll++)
					{

						rkml = sqrt(pow(a_arr[kk], 2.0) + pow(bo_arr[mm], 2.0) + pow(co_arr[ll], 2.0));
						mult_fact = pow(-1, (kk + mm + ll));//no plus 1 because number from 0

						if ((std::abs(bo_arr[mm] + rkml) < eps2) || std::abs(rkml) < eps2)
						{
							term1 = 0;
						}
						else
						{
							term1_1 = 3.0 * pow(a_arr[kk], 2.0)*pow(co_arr[ll], 2.0) - pow(co_arr[ll], 4.0);
							term1_2 = 3.0 * (pow(a_arr[kk], 2.0) - pow(co_arr[ll], 2.0))*(pow(a_arr[kk], 2.0) + pow(bo_arr[mm], 2.0) + pow(co_arr[ll], 2.0) + bo_arr[mm] * rkml)*std::log(bo_arr[mm] + rkml);
							term1_3 = 6.0 * rkml*(bo_arr[mm] + rkml);
							term1 = (term1_1 + term1_2) / term1_3;
						}
						if (std::abs(rkml) < eps2)
						{
							term2 = 0;
						}
						else
						{
							term2 = bo_arr[mm] * (3.0 * pow(a_arr[kk], 2.0) - pow(bo_arr[mm], 2.0)) / (6.0 * rkml);
						}
						if ((std::abs(a_arr[kk] + rkml) < eps2) || std::abs(rkml) < eps2)
						{
							term3_1 = 0;
						}
						else
						{
							term3_1 = (a_arr[kk] * bo_arr[mm] * pow(co_arr[ll], 2.0)) / (rkml*(a_arr[kk] + rkml));
						}
						if (std::abs(a_arr[kk] + rkml) < eps2)
						{
							term3_2 = 0;
						}
						else
						{
							term3_2 = a_arr[kk] * bo_arr[mm] * std::log(a_arr[kk] + rkml);
						}
						term3 = term3_1 + term3_2;

						if (std::abs(rkml) < eps2)
						{
							term4 = 0;
						}
						else
						{
							term4 = bo_arr[mm] * (pow(a_arr[kk], 2.0) + pow(bo_arr[mm], 2.0) + 2 * pow(co_arr[ll], 2.0)) / (3.0*rkml);
						}

						if ((std::abs(rkml) < eps2) || (std::abs(pow(a_arr[kk], 2.0) + pow(co_arr[ll], 2.0)) < eps2))
						{
							term5 = 0;
						}
						else
						{
							term5 = (pow(a_arr[kk], 4.0)*bo_arr[mm]) / (6.0 * (pow(a_arr[kk], 2.0) + pow(co_arr[ll], 2.0))*rkml);
						}

						if ((std::abs(rkml) < eps2) || (std::abs(pow(bo_arr[mm], 2) + pow(co_arr[ll], 2)) < eps2))
						{
							term6 = 0;
						}
						else
						{
							term6 = (pow(a_arr[kk], 2.0)*pow(bo_arr[mm], 3.0)) / (2.0 * (pow(bo_arr[mm], 2.0) + pow(co_arr[ll], 2.0))*rkml);
						}

						term7_a = (a_arr[kk] * co_arr[ll] / 2.0);
						if ((std::abs(rkml) < eps2) || (std::abs(pow(bo_arr[mm], 2.0) + pow(co_arr[ll], 2.0)) < eps2) || (std::abs(pow(a_arr[kk], 2.0) + pow(co_arr[ll], 2.0)) < eps2))
						{
							term7_b = 0;
						}
						else
						{
							term7_b1 = (a_arr[kk] * bo_arr[mm] * co_arr[ll] * (pow(a_arr[kk], 2.0) + pow(bo_arr[mm], 2.0) + 2.0 * pow(co_arr[ll], 2.0)));
							term7_b2 = (pow(a_arr[kk], 2.0) + pow(co_arr[ll], 2.0))*(pow(bo_arr[mm], 2.0) + pow(co_arr[ll], 2.0))*rkml;
							term7_b = term7_b1 / term7_b2;
						}
						if (std::abs(rkml*co_arr[ll]) < eps2)
						{
							term7_c = 0;
						}
						else
						{
							term7_c = 2.0 * std::atan(a_arr[kk] * bo_arr[mm] / (co_arr[ll] * rkml));
						}
						term7 = term7_a * (-term7_b + term7_c);

						int_res += (mult_fact*(term1 + term2 + term3 - term4 - term5 - term6 - term7));
					}
				}
			}
			if (std::abs((unit_normal_obs[0] + unit_normal_obs[1] + unit_normal_obs[2]) - sign(CC_w_sgn)) > eps2)
			{
				int_res = 0 - int_res;
			}
		}
	}

	return int_res;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{

	double dx = mxGetScalar(prhs[0]);
	double* src_pnl_loc = (double*)mxGetPr(prhs[1]);
	double* obs_pnl_loc = (double*)mxGetPr(prhs[2]);
	double* src_pnl_normal = (double*)mxGetPr(prhs[3]);
	double* obs_pnl_normal = (double*)mxGetPr(prhs[4]);
	double* smpl_pnts = (double*)mxGetPr(prhs[5]);
	double* smpl_wghts = (double*)mxGetPr(prhs[6]);
	int len_smpl = mxGetScalar(prhs[7]);
	int num_diff_pnts = mxGetScalar(prhs[8]);
	
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* int_rel_temp = (double*)mxGetPr(plhs[0]);

	double src_cen_temp[3];
	double obs_cen_temp[3];
	double unit_normal_src_temp[3];
	double unit_normal_obs_temp[3];
	for (int ii = 0; ii < 3; ii++)
	{
		src_cen_temp[ii] = *(src_pnl_loc + ii);
		obs_cen_temp[ii] = *(obs_pnl_loc + ii);
		unit_normal_src_temp[ii] = *(src_pnl_normal + ii);
		unit_normal_obs_temp[ii] = *(obs_pnl_normal + ii);
	}

	double int_res = 0;
	double temp_pnl[3];
	temp_pnl[0] = obs_pnl_loc[0] - src_pnl_loc[0];
	temp_pnl[1] = obs_pnl_loc[1] - src_pnl_loc[1];
	temp_pnl[2] = obs_pnl_loc[2] - src_pnl_loc[2];
	if (norm(temp_pnl) > 6.5 * dx)
	{
		int_res = finite_difference_coefficients(dx, num_diff_pnts, src_pnl_loc, obs_pnl_loc, src_pnl_normal, obs_pnl_normal, smpl_wghts, smpl_pnts, len_smpl);
	}
	else
	{
		int_res = analytical_derivation_1overR( dx, src_pnl_loc, obs_pnl_loc, src_pnl_normal, obs_pnl_normal);
	}
	int_rel_temp[0] = int_res;
}


