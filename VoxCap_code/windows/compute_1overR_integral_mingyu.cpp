
#include <string.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <malloc.h>
#include "compute_1overR_integral_mingyu.h"
//#include "compute_1overR_quad_mingyu.h"
//#include "compute_1overR_analy_mingyu.h"
#include "mex.h"
using namespace std;
template<class T>
double norm(T& arr)
{
	return sqrt(pow(arr[0], 2) + pow(arr[1], 2) + pow(arr[2], 2));
}

double compute_1overR_quad_mingyu(double src_cen[3], double obs_cen[3], double unit_normal_src[3], double unit_normal_obs[3], double smpl_wghts[], double smpl_pnts_temp[], int len_smpl)
{

	double arr_smpl[4096 * 4];
	double* smpl_pnts = &arr_smpl[0];

	for (int ii = 0; ii < len_smpl * 4; ii++)
	{
		arr_smpl[ii] = *(smpl_pnts_temp + ii);
	}

	int plane_src = 0; int plane_obs = 0; int form_type = 0; int ii;
	double y = 0.0;
	double int_res = 0.0;

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



double compute_1overR_analy_mingyu(double dx, double src_cen[3], double obs_cen[3], double unit_normal_src[3], double unit_normal_obs[3])
{

	double int_res = 0.0;
	int fl_self_term; int kk; int mm; int ll;
	double arg_log;
	double term1; double term2; double term3; double term4; double term5; double term6; double term7; double y;
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
		//arg_log = (dx + sqrt(2 * pow(dx, 2))) / dx;
		arg_log = 1 + sqrt(2);
		term1 = 3 * pow(dx, 3)*std::log(arg_log);
		term2 = pow(2, 1.5)*pow(dx, 3);
		//term2 = pow((2 * pow(dx, 2)), 3 / 2);
		term3 = term1;
		term4 = 2 * pow(dx, 3);
		int_res = (term1 - term2 + term3 + term4) * 2 / 3;
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
			}

			else if (plane_src == 2)
			{
				z = std::abs(obs_cen[1] - src_cen[1]);
				aij = std::abs(obs_cen[2] - src_cen[2]);
				bij = std::abs(obs_cen[0] - src_cen[0]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bp_arr[0] = bij - dx; bp_arr[1] = bij; bp_arr[2] = bij + dx; bp_arr[3] = bij;
			}
			else if (plane_src == 3)
			{
				z = std::abs(obs_cen[0] - src_cen[0]);
				aij = std::abs(obs_cen[1] - src_cen[1]);
				bij = std::abs(obs_cen[2] - src_cen[2]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bp_arr[0] = bij - dx; bp_arr[1] = bij; bp_arr[2] = bij + dx; bp_arr[3] = bij;
			}

			double a2minusz2over2[4]; double b2minz2over2[4]; double b2minus2z2[4]; double bz[4];
            a2minusz2over2[0]= (pow(a_arr[0], 2) - pow(z, 2)) / 2.0;   a2minusz2over2[1]= (pow(a_arr[1], 2) - pow(z, 2)) / 2.0;
			a2minusz2over2[2]= (pow(a_arr[2], 2) - pow(z, 2)) / 2.0;   a2minusz2over2[3]= (pow(a_arr[3], 2) - pow(z, 2)) / 2.0; 
			b2minz2over2[0]=(pow(bp_arr[0], 2) - pow(z, 2)) / 2.0; b2minz2over2[1]=(pow(bp_arr[1], 2) - pow(z, 2)) / 2.0;
			b2minz2over2[2]=(pow(bp_arr[2], 2) - pow(z, 2)) / 2.0; b2minz2over2[3]=(pow(bp_arr[3], 2) - pow(z, 2)) / 2.0;

			b2minus2z2[0] = pow(bp_arr[0], 2) - 2.0 * pow(z, 2); b2minus2z2[1]= pow(bp_arr[1], 2) - 2.0 * pow(z, 2);
			b2minus2z2[2] = pow(bp_arr[2], 2) - 2.0 * pow(z, 2); b2minus2z2[3]= pow(bp_arr[3], 2) - 2.0 * pow(z, 2);

			bz[0] =bp_arr[0] * z;  bz[1] =bp_arr[1] * z; bz[2] = bp_arr[2] * z; bz[3] =bp_arr[3] * z;

			double rkm; double mult_fact;
			for (kk = 0; kk < 4; kk++)
			{
				for (mm = 0; mm < 4; mm++)
				{
					rkm = sqrt(pow(a_arr[kk], 2) + pow(bp_arr[mm], 2) + pow(z, 2));
					mult_fact = pow(-1, (mm + kk));
					if (std::abs(a_arr[kk] + rkm) < 1e-37)
					{
						term1 = 0;
					}
					else
					{
						term1 = b2minz2over2[mm] * a_arr[kk] * std::log(a_arr[kk] + rkm);
					}

					if (std::abs(bp_arr[mm] + rkm) < 1e-37)
					{
						term2 = 0;
					}
					else
					{
						term2 = a2minusz2over2[kk] * bp_arr[mm] * std::log(bp_arr[mm] + rkm);
					}
					term3 = (1.0 / 6.0)*(b2minus2z2[mm] + pow(a_arr[kk], 2))*rkm;
					if (std::abs(z) < 1e-37)
					{
						term4 = 0;
					}
					else
					{
						term4 = bz[mm] * a_arr[kk] * std::atan((a_arr[kk] * bp_arr[mm]) / (rkm*z));
					}
					int_res += (mult_fact * (term1 + term2 - term3 - term4));
				}
			}

		}
		else if (form_type == 2)
		{
			if ((plane_src == 1 && plane_obs == 2) || (plane_src == 2 && plane_obs == 1))
			{
				aij = std::abs(obs_cen[0] - src_cen[0]);
				bij = std::abs(obs_cen[1] - src_cen[1]);
				cij = std::abs(obs_cen[2] - src_cen[2]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
			}
			else if ((plane_src == 1 && plane_obs == 3) || (plane_src == 3 && plane_obs == 1))
			{
				aij = std::abs(obs_cen[1] - src_cen[1]);
				bij = std::abs(obs_cen[0] - src_cen[0]);
				cij = std::abs(obs_cen[2] - src_cen[2]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
			}
			else if ((plane_src == 2 && plane_obs == 3) || (plane_src == 3 && plane_obs == 2))
			{
				aij = std::abs(obs_cen[2] - src_cen[2]);
				bij = std::abs(obs_cen[0] - src_cen[0]);
				cij = std::abs(obs_cen[1] - src_cen[1]);
				a_arr[0] = aij - dx; a_arr[1] = aij; a_arr[2] = aij + dx; a_arr[3] = aij;
				bo_arr[0] = bij + dx / 2.0; bo_arr[1] = bij - dx / 2.0;
				co_arr[0] = cij + dx / 2.0; co_arr[1] = cij - dx / 2.0;
			}

			double a2over2[4] = { pow(a_arr[0], 2.0) / 2.0, pow(a_arr[1], 2.0) / 2.0, pow(a_arr[2], 2.0) / 2.0, pow(a_arr[3], 2.0) / 2.0 };
			double c2over6[2] = { pow(co_arr[0], 2.0) / 6.0, pow(co_arr[1], 2.0) / 6.0 };
			double b2over6[2] = { pow(bo_arr[0], 2.0) / 6.0, pow(bo_arr[1], 2.0) / 6.0 };
			double a3over6[4] = { pow(a_arr[0], 3.0) / 6.0, pow(a_arr[1], 3.0) / 6.0, pow(a_arr[2], 3.0) / 6.0, pow(a_arr[3], 3.0) / 6.0 };
			double b2over2[2] = { pow(bo_arr[0], 2.0) / 2.0, pow(bo_arr[1], 2.0) / 2.0 };
			double c2over2[2] = { pow(co_arr[0], 2.0) / 2.0, pow(co_arr[1], 2.0) / 2.0 };

			for (kk = 0; kk < 4; kk++)
			{
				for (mm = 0; mm < 2; mm++)
				{
					for (ll = 0; ll < 2; ll++)
					{
						double rkml; double mult_fact;
						rkml = sqrt(pow(a_arr[kk], 2.0) + pow(bo_arr[mm], 2.0) + pow(co_arr[ll], 2.0));
						mult_fact = pow(-1, (kk + mm + ll));
						if (std::abs(bo_arr[mm] + rkml) < 1e-37)
						{
							term1 = 0;
						}
						else
						{
							term1 = (a2over2[kk] - c2over6[ll])*co_arr[ll] * std::log(bo_arr[mm] + rkml);
						}
						if (std::abs(co_arr[ll] + rkml) < 1e-37)
						{
							term2 = 0;
						}
						else
						{
							term2 = (a2over2[kk] - b2over6[mm])*bo_arr[mm] * std::log(co_arr[ll] + rkml);
						}
						if (std::abs(a_arr[kk] + rkml) < 1e-37)
						{
							term3 = 0;
						}
						else
						{
							term3 = a_arr[kk] * bo_arr[mm] * co_arr[ll] * std::log(a_arr[kk] + rkml);
						}
						term4 = bo_arr[mm] * co_arr[ll] * rkml / 3.0;
						if (std::abs(a_arr[kk] * rkml) < 1e-37)
						{
							term5 = 0;
						}
						else
						{
							term5 = a3over6[kk] * std::atan((bo_arr[mm] * co_arr[ll]) / (rkml*a_arr[kk]));
						}
						if (std::abs(bo_arr[mm] * rkml) < 1e-37)
						{
							term6 = 0;
						}
						else
						{
							term6 = b2over2[mm] * a_arr[kk] * std::atan((a_arr[kk] * co_arr[ll]) / (rkml*bo_arr[mm]));
						}
						if (std::abs(co_arr[ll] * rkml) < 1e-37)
						{
							term7 = 0;
						}
						else
						{
							term7 = c2over2[ll] * a_arr[kk] * std::atan((a_arr[kk] * bo_arr[mm]) / (rkml*co_arr[ll]));
						}
						int_res += mult_fact * (term1 + term2 + term3 - term4 - term5 - term6 - term7);

					}
				}
			}

		}
	}
	return int_res;

}



//double compute_1overR_integral_mingyu(double &dx, double (&src_cen)[3],
//	double (&obs_cen)[3], double (&unit_normal_src)[3], double
//	(&unit_normal_obs)[3], double smpl_pnts_temp[], double smpl_wghts[], int &len_smpl)


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{

	double dx = mxGetScalar(prhs[0]);
	double* src_cen =(double*)mxGetPr(prhs[1]);
	double* obs_cen = (double*)mxGetPr(prhs[2]);
	double* unit_normal_src = (double*)mxGetPr(prhs[3]);
	double* unit_normal_obs = (double*)mxGetPr(prhs[4]);
	double* smpl_pnts_temp = (double*)mxGetPr(prhs[5]);
	double* smpl_wghts = (double*)mxGetPr(prhs[6]);
	int len_smpl = mxGetScalar(prhs[7]);
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double* int_rel_temp = (double*)mxGetPr(plhs[0]);






	double int_res = 0;
	double temp_pnl[3];
	temp_pnl[0] = obs_cen[0] - src_cen[0];
	temp_pnl[1] = obs_cen[1] - src_cen[1];
	temp_pnl[2] = obs_cen[2] - src_cen[2];
	if (norm(temp_pnl) > 100 * dx)
	{
		int_res = compute_1overR_quad_mingyu( src_cen, obs_cen,  unit_normal_src,  unit_normal_obs,  smpl_wghts,  smpl_pnts_temp,  len_smpl);
	}
	else
	{
		int_res = compute_1overR_analy_mingyu(dx,  src_cen,  obs_cen,  unit_normal_src,  unit_normal_obs);
	}
	//return int_res;
	int_rel_temp[0] = int_res;
}


