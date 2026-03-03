#ifndef test_3D_SquareToCircle_logcost_MonteCarlo
#define test_3D_SquareToCircle_logcost_MonteCarlo

#include <string>

#include "../QuasiMonteCarlo/Generic_3D_logcost_MonteCarlo.h"
#include "../SmallGrid/SmallSinkhorn_3D_MC.h"


string testname="3D_SquareToCircle_logcost_MonteCarlo";

// Spherical box bounds (in radians, set at runtime via main())
double src_theta_min, src_theta_max, src_phi_min, src_phi_max;
double tgt_theta_min, tgt_theta_max, tgt_phi_min, tgt_phi_max;


double P(double x[])
{
	double theta_c     = 0.5 * (src_theta_min + src_theta_max);
	double phi_c       = 0.5 * (src_phi_min   + src_phi_max);
	double sigma_theta = 0.5 * (src_theta_max  - src_theta_min);
	double sigma_phi   = 0.5 * (src_phi_max    - src_phi_min);

	double theta = acos(x[2]);
	double phi   = atan2(x[1], x[0]);
	if (phi < 0) phi += 2*PI;

	double dtheta = theta - theta_c;
	double dphi   = phi   - phi_c;
	if (dphi >  PI) dphi -= 2*PI;
	if (dphi < -PI) dphi += 2*PI;

	return exp(-0.5 * (dtheta*dtheta/(sigma_theta*sigma_theta)
	                 + dphi  *dphi  /(sigma_phi  *sigma_phi  )));
}


double Q(double y[])
{
	double theta_c     = 0.5 * (tgt_theta_min + tgt_theta_max);
	double phi_c       = 0.5 * (tgt_phi_min   + tgt_phi_max);
	double sigma_theta = 0.5 * (tgt_theta_max  - tgt_theta_min);
	double sigma_phi   = 0.5 * (tgt_phi_max    - tgt_phi_min);

	double theta = acos(y[2]);
	double phi   = atan2(y[1], y[0]);
	if (phi < 0) phi += 2*PI;

	double dtheta = theta - theta_c;
	double dphi   = phi   - phi_c;
	if (dphi >  PI) dphi -= 2*PI;
	if (dphi < -PI) dphi += 2*PI;

	return exp(-0.5 * (dtheta*dtheta/(sigma_theta*sigma_theta)
	                 + dphi  *dphi  /(sigma_phi  *sigma_phi  )));
}


#endif