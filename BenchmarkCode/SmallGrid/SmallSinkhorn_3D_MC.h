#ifndef SmallSinkhorn_3D_MC
#define SmallSinkhorn_3D_MC

#include "../QuasiMonteCarlo/Generic_3D_logcost_MonteCarlo.h"


double p_small[NK_small];
double q_small[NK_small];

double F_small[NK_small];
double G_small[NK_small];

double f_small[NK_small];
double g_small[NK_small];

double temp_small[NK_small];


void smallsinkhorn(int regvariable)
{

	fill_n(F_small,NK_small,1);
	fill_n(G_small,NK_small,1);

	fill_n(f_small,NK_small,0);
	fill_n(g_small,NK_small,0);

	double sum=0;
	for(int i=0; i<NK_small; i++)
	{
		p_small[i]=P(x_small[i]);
		sum+=p_small[i];
	}	

	for(int i=0; i<NK_small; i++)
	{
		p_small[i]/=sum;
	}


	sum=0;
	for(int i=0; i<NK_small; i++)
	{
		q_small[i]=Q(y_small[i]);
		sum+=q_small[i];
	}


	for(int i=0; i<NK_small; i++)
	{
		q_small[i]/=sum;
	}

	sum=0;
	for(int i=0; i<NK_small; i++)
	{
		for(int j=0; j<NK_small; j++)
		{
			temp_small[i]+=exp(-regvariable*Cost_Func(x_small[i],y_small[j]))*p_small[j];
		}
		sum+=temp_small[i];
	}

	for(int i=0; i<NK_small; i++)
	{
		temp_small[i]/=sum;
	}

	swap(temp_small,p_small);


	sum=0;
	for(int i=0; i<NK_small; i++)
	{
		for(int j=0; j<NK_small; j++)
		{
			temp_small[i]+=exp(-regvariable*Cost_Func(x_small[i],y_small[j]))*q_small[j];
		}
		sum+=temp_small[i];
	}

	for(int i=0; i<NK_small; i++)
	{
		temp_small[i]/=sum;
	}

	swap(temp_small,q_small);


	for(int k=1; k<regvariable; k+=pow(regvariable, 1./3))
	
	{

		for(int i=0; i<NK_small; i++)
		{
			if(p_small[i]==0) continue;
			double sum=0;
			double temp=0;

			for(int j=0; j<NK_small; j++)
			{
				temp=-k*(Cost_Func(x_small[i],y_small[j])-f_small[i]-g_small[j]);
				temp=exp(temp);
				sum+=temp*G_small[j]*q_small[j];

			}

			if(sum>0)F_small[i]=1./sum;
			else F_small[i]=1;

		
		}


		for(int j=0; j<NK_small; j++)
		{
			if (q_small[j]==0) continue;
			double sum=0;
			double temp=0;

			for(int i=0; i<NK_small; i++)
			{
				temp=-k*(Cost_Func(x_small[i],y_small[j])-f_small[i]-g_small[j]);
				temp=exp(temp);
				sum+=temp*F_small[i]*p_small[i];
			}

			if(sum>0)G_small[j]=1./sum;
			else G_small[j]=1;


		}

		for(int i=0; i<NK_small; i++)
		{
			f_small[i]+=log(F_small[i])/k;
			g_small[i]+=log(G_small[i])/k;

			F_small[i]=1;
			G_small[i]=1;


		}

	}
	
}

#endif