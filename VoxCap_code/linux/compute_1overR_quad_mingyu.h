#pragma once
#ifndef COMPUTE_1OVERR_QUAD_MINGYU_H
#define COMPUTE_1OVERR_QUAD_MINGYU_H
#include "mex.h"
extern double compute_1overR_quad_mingyu(double src_cen[3], double obs_cen[3], double unit_normal_src[3], double unit_normal_obs[3], double smpl_wghts[], double smpl_pnts_temp[], int len_smpl);

#endif
