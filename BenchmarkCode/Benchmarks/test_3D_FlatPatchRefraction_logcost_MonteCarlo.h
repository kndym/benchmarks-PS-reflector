#ifndef test_3D_FlatPatchRefraction_logcost_MonteCarlo
#define test_3D_FlatPatchRefraction_logcost_MonteCarlo

#include <string>
#include <random>
#include <cmath>

#include "../QuasiMonteCarlo/Generic_3D_refractioncost_MonteCarlo.h"
#include "../SmallGrid/SmallSinkhorn_3D_MC.h"


string testname = "3D_FlatPatchRefraction_logcost_MonteCarlo";

// Refractor problem from paper (kappa = 0.6):
//   Source  Omega  = { (1,phi,theta) : pi/12 < phi < pi/4,  pi/12 < theta < pi/3  }
//   Target  Omega* = { (1,phi,theta) : pi/10 < phi < pi/5,  pi/10 < theta < pi/5  }
// Both patches are on the upper hemisphere (theta < pi/2).
// Densities are flat (indicator = 1 inside patch, 0 outside).
//
// Notebook Step 1 defaults (values are multiples of pi):
//   SRC_THETA_MIN=1/12  SRC_THETA_MAX=1/3   SRC_PHI_MIN=1/12  SRC_PHI_MAX=1/4
//   TGT_THETA_MIN=1/10  TGT_THETA_MAX=1/5   TGT_PHI_MIN=1/10  TGT_PHI_MAX=1/5
//
// IMPORTANT: in the notebook Step 3, change y_pts to upper=True.

// Runtime patch bounds — set by main() from argv (same as other benchmarks)
double src_theta_min = PI / 12, src_theta_max = PI / 3;
double src_phi_min   = PI / 12, src_phi_max   = PI / 4;
double tgt_theta_min = PI / 10, tgt_theta_max = PI / 5;
double tgt_phi_min   = PI / 10, tgt_phi_max   = PI / 5;


// Fills x[NK] with source-patch points and y[NK] with target-patch points,
// both uniformly distributed on the sphere surface within their respective patches.
// Sampling: phi ~ Uniform[phi_min, phi_max], cos(theta) ~ Uniform[cos(theta_max), cos(theta_min)]
void generate_patch_points()
{
	std::mt19937 rng(42);
	std::uniform_real_distribution<double> U(0.0, 1.0);

	double cos_src_lo = cos(src_theta_max);
	double cos_src_hi = cos(src_theta_min);
	for (int i = 0; i < NK; i++) {
		double u   = cos_src_lo + (cos_src_hi - cos_src_lo) * U(rng);
		double phi = src_phi_min + (src_phi_max - src_phi_min) * U(rng);
		double s   = sqrt(1.0 - u * u);
		x[i][0] = s * cos(phi);
		x[i][1] = s * sin(phi);
		x[i][2] = u;
	}

	double cos_tgt_lo = cos(tgt_theta_max);
	double cos_tgt_hi = cos(tgt_theta_min);
	for (int i = 0; i < NK; i++) {
		double u   = cos_tgt_lo + (cos_tgt_hi - cos_tgt_lo) * U(rng);
		double phi = tgt_phi_min + (tgt_phi_max - tgt_phi_min) * U(rng);
		double s   = sqrt(1.0 - u * u);
		y[i][0] = s * cos(phi);
		y[i][1] = s * sin(phi);
		y[i][2] = u;
	}

	// The shared point-cloud header's small arrays are otherwise populated for
	// the square/circle reflector benchmark.  Refill them for this patch so
	// SmallSinkhorn evaluates the same P/Q supports as the main grid.
	double cos_src_small_lo = cos_src_lo;
	double cos_src_small_hi = cos_src_hi;
	for (int i = 0; i < NK_small; i++) {
		double u   = cos_src_small_lo + (cos_src_small_hi - cos_src_small_lo) * U(rng);
		double phi = src_phi_min + (src_phi_max - src_phi_min) * U(rng);
		double s   = sqrt(1.0 - u * u);
		x_small[i][0] = s * cos(phi);
		x_small[i][1] = s * sin(phi);
		x_small[i][2] = u;
	}

	double cos_tgt_small_lo = cos_tgt_lo;
	double cos_tgt_small_hi = cos_tgt_hi;
	for (int i = 0; i < NK_small; i++) {
		double u   = cos_tgt_small_lo + (cos_tgt_small_hi - cos_tgt_small_lo) * U(rng);
		double phi = tgt_phi_min + (tgt_phi_max - tgt_phi_min) * U(rng);
		double s   = sqrt(1.0 - u * u);
		y_small[i][0] = s * cos(phi);
		y_small[i][1] = s * sin(phi);
		y_small[i][2] = u;
	}
}


double P(double x[])
{
	double theta = acos(x[2]);
	double phi   = atan2(x[1], x[0]);
	if (phi < 0) phi += 2 * PI;

	if (theta >= src_theta_min && theta <= src_theta_max &&
	    phi   >= src_phi_min   && phi   <= src_phi_max)
		return 1.0;

	return 0.0;
}


double Q(double y[])
{
	double theta = acos(y[2]);
	double phi   = atan2(y[1], y[0]);
	if (phi < 0) phi += 2 * PI;

	if (theta >= tgt_theta_min && theta <= tgt_theta_max &&
	    phi   >= tgt_phi_min   && phi   <= tgt_phi_max)
		return 1.0;

	return 0.0;
}


#endif
