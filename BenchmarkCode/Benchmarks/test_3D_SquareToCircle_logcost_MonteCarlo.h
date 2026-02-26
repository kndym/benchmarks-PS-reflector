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
	double theta = acos(x[2]);           // polar angle [0, π]
	double phi   = atan2(x[1], x[0]);   // azimuthal [-π, π]
	if(phi < 0) phi += 2*PI;            // normalize to [0, 2π]

	if(theta >= src_theta_min && theta <= src_theta_max &&
	   phi   >= src_phi_min   && phi   <= src_phi_max)
		return 1;
	return 0;
}


double Q(double y[])
{
	double theta = acos(y[2]);
	double phi   = atan2(y[1], y[0]);
	if(phi < 0) phi += 2*PI;

	if(theta >= tgt_theta_min && theta <= tgt_theta_max &&
	   phi   >= tgt_phi_min   && phi   <= tgt_phi_max)
		return 1;
	return 0;
}


#endif