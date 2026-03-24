#ifndef test_3D_SquareToCircle_logcost_MonteCarlo
#define test_3D_SquareToCircle_logcost_MonteCarlo

#include <string>

#include "../QuasiMonteCarlo/Generic_3D_logcost_MonteCarlo.h"
#include "../SmallGrid/SmallSinkhorn_3D_MC.h"


string testname="3D_SquareToCircle_logcost_MonteCarlo";




double P(double x[])
{
	double proj1=x[0]/(1+x[2]);
	double proj2=x[1]/(1+x[2]);


	if(-0.5<proj1 && proj1<0.5 && -0.5<proj2 && proj2<0.5)
	{
		return 1;	
	}

	return 0;
}


double Q(double y[])
{
	double proj1=y[0]/(1-y[2]);
	double proj2=y[1]/(1-y[2]);


	if(proj1*proj1+proj2*proj2<=0.25)
	{
		return 1;
	}

	return 0;
}


#endif