#ifndef Generic_3D_logcost_MonteCarlo
#define Generic_3D_logcost_MonteCarlo

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cfloat>
#include <math.h>
#include <vector>

#include "../QuasiMonteCarlo/MonteCarlo_Pointcloud_3D_128.h"

const double PI = 3.14159265358979323846;

using namespace std;

	//according to dimension, the way x and y are processed will change.

int SpatialResolution;


double p[NK];
double q[NK];
double Original_p[NK];
double Original_q[NK];
double logp[NK];
double logq[NK];


//I'm putting here 128 instead of sqrt(NK) as this is greatest number and I want all of them to have uniform setting
const int MeshGridResolution=128+1;
const int DestinationResolution=65;
const int FinalGridResolution=2*128+1;
const int FinalGrid=FinalGridResolution*FinalGridResolution;
double x_regular[FinalGrid][dim];
double f_regular[FinalGrid];
double Ref_regular[FinalGrid][dim];

double Regular_side[FinalGridResolution];
double Y_pushed_projected[NK][dim];
double Y_Traced[DestinationResolution][DestinationResolution];

 //p,q measures, x,y mesh points,

double norm1=0,norm2=0;
int ind1=0,ind2=0;



vector < vector <double> > u(NK, vector<double> (NK,0));
vector < vector <double> > xproj (NK, vector<double> (dim-1,0));

double P(double x[]);
double Q(double x[]);



int getk(int Nk)
{
	return sqrt(Nk);
}


double Cost_Func(double x[], double y[])
{
	double temp=1;
	for(int d=0; d<dim; d++)
	{
		temp-=x[d]*y[d];
	}
	return -log(temp); //cost from WANG paper
}


double dist(double x[], double y[])
{
	double temp=0;
	for(int d=0; d<dim; d++)
	{
		temp+=(x[d]*y[d]);
	}
	return acos(temp); //Spherical distance
}


void discretization()
{
	fill_n(Original_p,NK,0);
	fill_n(Original_q,NK,0);
	SpatialResolution=getk(NK);
	//square to square, but if densities change I can give any shape I want


	for(int i=0; i<FinalGridResolution; i++)
	{
		for(int j=0; j<FinalGridResolution; j++)
		{
			double X=(-0.6+1.2*(double)i/(double)(FinalGridResolution-1));
			double Y=(-0.6+1.2*(double)j/(double)(FinalGridResolution-1));
			double N2=X*X+Y*Y;

			if(i==0)
			{
				Regular_side[j]=Y;
			}
			x_regular[i*FinalGridResolution+j][0]=2*X/(1+N2);
			x_regular[i*FinalGridResolution+j][1]=2*Y/(1+N2);
			x_regular[i*FinalGridResolution+j][2]=(1-N2)/(1+N2);		
		}
	}
}

void fill_p()
{
	int count=0;
	double sum=0;
	//filling uniform distribution on square x
	for(int i=0; i<NK; i++)
	{
		p[i]=P(x[i]);
		sum+=p[i];

		if(p[i]>0)
		{			
			count++;
		}
	}


	for(int i=0; i<NK; i++)
	{
		p[i]/=sum;
		logp[i]=log(p[i]);
	}
	cout<<"Number of Points in Source:		"<<count<<endl;

}

void fill_q()
{
	//circular density on the square y
	int count=0;
	int count1=0;
	double sum=0;
	for(int i=0; i<NK; i++)
	{
		q[i]=Q(y[i]);
		sum+=q[i];

		if(q[i]>0)
		{
			count++;			
		}
	}

	for(int i=0; i<NK; i++)
	{
		q[i]/=sum;
		logq[i]=log(q[i]);
	}
	cout<<"Number of Points in Destination:		"<<count<<endl;
}






double L2_norm(double x[])
{
	double sum=0;
	for(int i=0; i<dim; i++)
	{
		sum+=x[i]*x[i];
	}
	return sqrt(sum);
}


void MeshGridDestinationDensity(string outputname)
{
	string command=outputname+"/Y_MeshGrid.txt";
	freopen(command.c_str(),"w",stdout);


	for (int i=0; i<MeshGridResolution; i++)
	{

		for(int j=0; j<MeshGridResolution; j++)
		{
			double X=(-0.6+1.2*i/(MeshGridResolution-1));
			double Y=(-0.6+1.2*j/(MeshGridResolution-1));
			double N2=X*X+Y*Y;

			double vec[3];
			vec[0]=2*X/(1+N2);
			vec[1]=2*Y/(1+N2);
			vec[2]=-1*(1-N2)/(1+N2);
			
			printf("%.*e ",	DECIMAL_DIG, Q(vec)*4/((1+X*X+Y*Y)*(1+X*X+Y*Y)));

		}	
		printf("\n");	
	}

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);

	command=outputname+"/Y_DiscreteMesh.txt";
	freopen(command.c_str(),"w",stdout);


	for (int i=0; i<MeshGridResolution; i++)
	{

		for(int j=0; j<MeshGridResolution; j++)
		{
			double X=(-0.6+1.2*i/(MeshGridResolution-1));
			double Y=(-0.6+1.2*j/(MeshGridResolution-1));
			double N2=X*X+Y*Y;

			double vec[3];
			vec[0]=2*X/(1+N2);
			vec[1]=2*Y/(1+N2);
			vec[2]=-1*(1-N2)/(1+N2);
			if(Q(vec)==0) continue;
			printf("%.*e %.*e %.*e \n",	DECIMAL_DIG, X, DECIMAL_DIG, Y, DECIMAL_DIG, Q(vec)*4/((1+X*X+Y*Y)*(1+X*X+Y*Y)));
		}	
	}

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);
}


void MeshGridSourceDensity(string outputname)
{
	string command=outputname+"/X_MeshGrid.txt";
	freopen(command.c_str(),"w",stdout);


	for (int i=0; i<MeshGridResolution; i++)
	{

		for(int j=0; j<MeshGridResolution; j++)
		{
			double X=(-0.6+1.2*i/(MeshGridResolution-1));
			double Y=(-0.6+1.2*j/(MeshGridResolution-1));
			double N2=X*X+Y*Y;

			double vec[3];
			vec[0]=2*X/(1+N2);
			vec[1]=2*Y/(1+N2);
			vec[2]=1*(1-N2)/(1+N2);
			
			printf("%.*e ",	DECIMAL_DIG, P(vec)*4/((1+X*X+Y*Y)*(1+X*X+Y*Y)));

		}	
		printf("\n");	
	}

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);

	command=outputname+"/X_DiscreteMesh.txt";
	freopen(command.c_str(),"w",stdout);


	for (int i=0; i<MeshGridResolution; i++)
	{

		for(int j=0; j<MeshGridResolution; j++)
		{
			double X=(-0.6+1.2*i/(MeshGridResolution-1));
			double Y=(-0.6+1.2*j/(MeshGridResolution-1));
			double N2=X*X+Y*Y;

			double vec[3];
			vec[0]=2*X/(1+N2);
			vec[1]=2*Y/(1+N2);
			vec[2]=1*(1-N2)/(1+N2);
			
			if(P(vec)==0) continue;
			printf("%.*e %.*e %.*e \n",	DECIMAL_DIG, X, DECIMAL_DIG, Y, DECIMAL_DIG, P(vec)*4/((1+X*X+Y*Y)*(1+X*X+Y*Y)));
		}	
		printf("\n");	
	}

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);
}

#endif