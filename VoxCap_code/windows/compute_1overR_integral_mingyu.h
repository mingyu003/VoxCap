
#ifndef COMPUTE_1OVERR_INTEGRAL_MINGYU_H
#define COMPUTE_1OVERR_INTEGRAL_MINGYU_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "mex.h"
extern double compute_1overR_integral_mingyu(double &dx, double (&src_pnl_loc)[3],
	 double (&obs_pnl_loc)[3], double (&src_pnl_normal)[3], double
	(&obs_pnl_normal)[3], double smpl_pnts[], double smpl_wghts[], int &len_smpl);

#endif

#pragma once
