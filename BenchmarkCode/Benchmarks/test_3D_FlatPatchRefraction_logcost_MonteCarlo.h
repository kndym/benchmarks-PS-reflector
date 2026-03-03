#ifndef test_3D_FlatPatchRefraction_logcost_MonteCarlo
#define test_3D_FlatPatchRefraction_logcost_MonteCarlo

#include <string>

#include "../QuasiMonteCarlo/Generic_3D_refractioncost_MonteCarlo.h"
#include "../SmallGrid/SmallSinkhorn_3D_MC.h"


string testname = "3D_FlatPatchRefraction_logcost_MonteCarlo";

// Refractor problem from paper (kappa = 0.6):
//   Source  Omega  = { (1,phi,theta) : pi/12 < phi < pi/4,  pi/12 < theta < pi/3  }
//   Target  Omega* = { (1,phi,theta) : pi/10 < phi < pi/5,  pi/10 < theta < pi/5  }
// Both patches are on the upper hemisphere (theta < pi/2).
// Densities are flat (indicator = 1 inside patch, 0 outside).
//
// IMPORTANT: regenerate the y point cloud with upper=True in the notebook
// (change gen_pts(NK, upper=False) -> gen_pts(NK, upper=True) for y_pts).

// Source patch bounds
static const double SRC_THETA_MIN = PI / 12.0;   // pi/12
static const double SRC_THETA_MAX = PI / 3.0;    // pi/3
static const double SRC_PHI_MIN   = PI / 12.0;   // pi/12
static const double SRC_PHI_MAX   = PI / 4.0;    // pi/4

// Target patch bounds
static const double TGT_THETA_MIN = PI / 10.0;   // pi/10
static const double TGT_THETA_MAX = PI / 5.0;    // pi/5
static const double TGT_PHI_MIN   = PI / 10.0;   // pi/10
static const double TGT_PHI_MAX   = PI / 5.0;    // pi/5


double P(double x[])
{
	double theta = acos(x[2]);
	double phi   = atan2(x[1], x[0]);
	if (phi < 0) phi += 2 * PI;

	if (theta >= SRC_THETA_MIN && theta <= SRC_THETA_MAX &&
	    phi   >= SRC_PHI_MIN   && phi   <= SRC_PHI_MAX)
		return 1.0;

	return 0.0;
}


double Q(double y[])
{
	double theta = acos(y[2]);
	double phi   = atan2(y[1], y[0]);
	if (phi < 0) phi += 2 * PI;

	if (theta >= TGT_THETA_MIN && theta <= TGT_THETA_MAX &&
	    phi   >= TGT_PHI_MIN   && phi   <= TGT_PHI_MAX)
		return 1.0;

	return 0.0;
}


#endif
