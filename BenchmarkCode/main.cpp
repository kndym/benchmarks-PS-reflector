#include 	"Benchmarks/test_3D_SquareToCircle_logcost_MonteCarlo.h"			//thats where p,q,x,y and cost come from
//  Choice of the testcase, use the desired name of the testcase instead of **** in the #include "****"

//Benchmark testcases:

//			"Benchmarks/test_3D_SquareToCircle_logcost_MonteCarlo.h"
//			"Benchmarks/test_3D_SquareToTwoGaussSide_logcost_MonteCarlo.h"



#include <ctime>
#include <iomanip>


//	Functions for computing reflection from the reflector computed on the regular (rectangular instead of monte-carlo) grid

#include "PushForward/Pushforward_of_RefRegular.h"


// --- MKL replacement functions (standard C math, no MKL dependency) ---
inline void vdMul(int n, const double* a, const double* b, double* r) {
	for(int i = 0; i < n; i++) r[i] = a[i] * b[i];
}
inline void vdLn(int n, const double* a, double* r) {
	for(int i = 0; i < n; i++) r[i] = log(a[i]);
}
inline void vdExp(int n, const double* a, double* r) {
	for(int i = 0; i < n; i++) r[i] = exp(a[i]);
}
inline void vdAdd(int n, const double* a, const double* b, double* r) {
	for(int i = 0; i < n; i++) r[i] = a[i] + b[i];
}
inline int cblas_idamax(int n, const double* x, int /*inc*/) {
	int idx = 0;
	for(int i = 1; i < n; i++) if(fabs(x[i]) > fabs(x[idx])) idx = i;
	return idx;
}
inline void cblas_dscal(int n, double alpha, double* x, int /*inc*/) {
	for(int i = 0; i < n; i++) x[i] *= alpha;
}


string name, outputname;
string command;

//	Variable for keeping track of change of potential f between iterations
double maxdif=-1;


//	variables for kantorovich potentials in sinkhorn iterations
// 	f,g potentials, F,G exponents of potentials 

double F[NK];
double G[NK];

double f[NK];
double g[NK];

//	Same as above, for diagonal sinkhorn iterations

double F_id[NK];
double G_id[NK];

double f_id[NK];
double g_id[NK];


//	Variables for builging reflector

double fc[NK];
double gc[NK];

double R[NK];
double Ref[NK][dim];
double Refc[NK][dim];


//	Temprorary containers for swapping intermediate values during iterations

double K_oneline[NK];

double tempvecNK_1[NK];
double tempvecNK_2[NK];
double tempvecNK_3[NK];



//	Treshold, below which the values will not be summed and the index will be declared as "useless" using index bounds

double treshold=log(1.e-6/NK);


//	Bounds for Indexes outside of which the summation should not happen in sinkhron iteration, as values are below above treshold

int IndJ_begin[NK];
int IndJ_end[NK];

int IndI_begin[NK];
int IndI_end[NK];

int IndId_begin[NK];
int IndId_end[NK];


//	Function for filling all declared variables with initial values

void fillall();


//	one sinkhorn iteration

void Sinkhorn_axb(int k, int iterator=1, double correctionterm =1.);


//	Absorbing values of potentials after every iteration, for stability of iterative steps
void absorbtion(int k);
void absorbtion_f_id(int k);
void absorbtion_g_id(int k);

//	One sinkhorn interation for "identity" computations

void Sinkhorn_identity_F_axb(int k, int iterator=1, double correctionterm =1.);
void Sinkhorn_identity_G_axb(int k, int iterator=1, double correctionterm =1.);


//	Sinkhorn algorithm without correcting entropic bias 
//and building reflector from the result
void do_sinkhorn_pure_axb(int k, double errortolerance);

//	(Sinkhorn divergence) Sinkhorn algorithm with correction of entropic bias by subtracting "identitiy" terms, 
//and building reflector from the result
void do_sinkhorn_subtracted_axb(int k, double errortolerance=1.E-12);


// Building reflector from the sinkhorn algorithm's outcome
void sink_to_Reflector_subtracted_axb(int k);


//	functions for c-transform computations for kantorovich potentials f and g

void GetCtransforms();
void Get_fc();
void Get_gc();



//	Printing functions saving varios information
void Printall_MY();

void print(double  vector [], int size=NK);
void print(int vector [], int size);
void print(double  vector [][dim], int size1=NK, int size2=dim);

//	Doing sinkhorn iterations on the small scale
void UseSmall(int regparam);


//	Printing the support of the desired target distribution.
void Get_Original_Y(string outputname);


//	Computing the total sinkhorn divergence cost of transport between source and target 
double Compute_TotalCost();



// Regularization parameter Reg_param=1/epsilon
double Reg_Param;




//	Value treshold cap and iteration number cap for sinkhorn iterations:
//	After reaching the final value of regularization parameter, iterations will continue with this final value
//till either maxdif<cap_treshold or number of iterations done at final scale reaches the cap_iterations

double cap_treshold=(1.e-5);
int cap_iteration=16;


//	Multiplier of spatial resolution k, to give the regularization parameter:
//			Reg_param = multiplier* Spatial resolution
int multiplier=8;

int main(int argc, char* argv[])
{
	if(argc != 9) {
		fprintf(stderr,
			"Usage: main src_th_min src_th_max src_ph_min src_ph_max"
			" tgt_th_min tgt_th_max tgt_ph_min tgt_ph_max\n"
			"(values as multiples of pi)\n");
		return 1;
	}
	src_theta_min = atof(argv[1])*PI;
	src_theta_max = atof(argv[2])*PI;
	src_phi_min   = atof(argv[3])*PI;
	src_phi_max   = atof(argv[4])*PI;
	tgt_theta_min = atof(argv[5])*PI;
	tgt_theta_max = atof(argv[6])*PI;
	tgt_phi_min   = atof(argv[7])*PI;
	tgt_phi_max   = atof(argv[8])*PI;

	clock_t Totaltime = clock();

	char name1[100];
	time_t rawtime;
  	struct tm * timeinfo;
  	time (&rawtime);
  	timeinfo = localtime (&rawtime);
  	strftime (name1,100,"_%Y_%m_%d_%H:%M:%S",timeinfo);

	name="Output"+string(name1)+"_"+to_string(dim)+"D_NK"+to_string(NK)+"_";
	outputname="Output"+string(name1);

	command="mkdir "+outputname;
	system(command.c_str());


	command=outputname+"/log.txt";
	freopen(command.c_str(),"w",stdout);

	string addname="_NumbOfIt__Steppow_1l3";




	fillall();
	Reg_Param=multiplier*getk(NK);



	cout<<"Dimension of matrix K: "<<NK<<endl;


	//do_sinkhorn_pure_axb(Reg_Param, 1./(Reg_Param));
	do_sinkhorn_subtracted_axb(Reg_Param, 1./(Reg_Param));
	

	GetCtransforms();


	Printall_MY();

	Get_Original_Y(outputname);

	MeshGridDestinationDensity(outputname);
	MeshGridSourceDensity(outputname);


	Pushforward_Ref_regular(outputname);



	Totaltime= clock()-Totaltime;
	printf("Total time in clicks %ld, time in seconds %f with change %.*e \n", Totaltime, (float)Totaltime/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);
	
	cout<<endl;

	command="Results_"+testname;
	system(("mkdir "+command).c_str());
	command="mv "+outputname+" "+command+"/"+name+addname;
	system(command.c_str());

	system("spd-say done ");


	return 0;
}


void fillall()
{
	//functions from test header
	discretization();
	fill_p();
	fill_q();

	fill_n(F,NK,1);
	fill_n(G,NK,1);

	fill_n(f,NK,0);
	fill_n(g,NK,0);

	fill_n(fc,NK,0);
	fill_n(gc,NK,0);

	fill_n(F_id,NK,1);
	fill_n(G_id,NK,1);

	fill_n(f_id,NK,0);
	fill_n(g_id,NK,0);



	fill_n(IndJ_begin,NK,0);
	fill_n(IndJ_end,NK,NK);
	fill_n(IndI_begin,NK,0);
	fill_n(IndI_end,NK, NK);

	fill_n(IndId_begin,NK,0);
	fill_n(IndId_end,NK, NK);


}



void Sinkhorn_identity_F_axb(int k, int iterator, double correctionterm)
{
	maxdif=-1;
	int maxind;

	vdMul(NK,F_id,p,tempvecNK_1);

	for(int i=0; i<NK; i++)
	{
		if(p[i]==0)continue;
		double sum=0;
		double temp=0;

		int ind_begin=0;
		int ind_end=0;

		bool first=true;

		for(int j=IndId_begin[i]; j<IndId_end[i]; j++)
		{
			temp=-k*(Cost_Func(x[i],y[j])-f_id[i]-f_id[j]);
			if(temp+logp[j]<treshold) continue;

			if(first)
			{
				ind_begin=j;
				first=false;
			}
			ind_end=j;

			temp=exp(temp);
			sum+=temp*tempvecNK_1[j];

		}
		IndId_begin[i]=ind_begin;
		IndId_end[i]=ind_end+1;

		if(sum>0)F_id[i]=1./sum;
		else F_id[i]=1;
	}

	//averaging previous value and new approximation (feyde's idea) seems to perform well in identity case
	
	vdLn(NK,F_id,tempvecNK_3);
	maxind=cblas_idamax(NK,tempvecNK_3,1);
	maxdif=fabs(tempvecNK_3[maxind]/k);

}



void Sinkhorn_identity_G_axb(int k, int iterator, double correctionterm)
{
	maxdif=-1;
	int maxind;

	vdMul(NK,G_id,q,tempvecNK_1);

	for(int j=0; j<NK; j++)
	{
		if(q[j]==0)continue;
		double sum=0;
		double temp=0;

		int ind_begin=0;
		int ind_end=0;

		bool first=true;

		for(int i=IndId_begin[j]; i<IndId_end[j]; i++)
		{
			temp=-k*(Cost_Func(x[i],y[j])-g_id[i]-g_id[j]);
			if(logq[i]+temp<treshold) continue;

			if(first)
			{
				ind_begin=i;
				first=false;
			}
			ind_end=i;

			temp=exp(temp);
			sum+=temp*tempvecNK_1[i];

		}
		IndId_begin[j]=ind_begin;
		IndId_end[j]=ind_end+1;

		if(sum>0)G_id[j]=1./sum;
		else G_id[j]=1;
	}
	//averaging previous value and new approximation (feyde's idea) seems to perform well in identity case
	vdLn(NK,G_id,tempvecNK_3);
	maxind=cblas_idamax(NK,tempvecNK_3,1);
	maxdif=fabs(tempvecNK_3[maxind]/k);

}



void Sinkhorn_axb(int k, int iterator, double correctionterm)
{
	int maxind=0;
	double maxval=0;

	
	vdMul(NK,F,p,tempvecNK_1);

	for(int j=0; j<NK; j++)
	{
		if(q[j]==0)continue;

		double sum=0;
		double temp=0;

		int ind_begin=0;
		int ind_end=0;

		bool first=true;

		for(int i=IndI_begin[j]; i<IndI_end[j]; i++)
		{
			temp=-k*(Cost_Func(x[i],y[j])-f[i]-g[j]);
			if(logp[i]+temp<treshold)continue;

			if(first)
			{
				ind_begin=i;
				first=false;
			}
			ind_end=i;

			temp=exp(temp);
			sum+=temp*tempvecNK_1[i];
		}

		IndI_begin[j]=ind_begin;
		IndI_end[j]=ind_end+1;

		if(sum>0)G[j]=1./sum;
		else G[j]=1;
	}

		vdMul(NK,G,q,tempvecNK_1);

	for(int i=0; i<NK; i++)
	{
		if(p[i]==0)continue;
		double sum=0;
		double temp=0;

		int ind_begin=0;
		int ind_end=0;

		bool first=true;

		for(int j=IndJ_begin[i]; j<IndJ_end[i]; j++)
		{
			temp=-k*(Cost_Func(x[i],y[j])-f[i]-g[j]);
			if(temp+logq[j]<treshold) continue;

			if(first)
			{
				ind_begin=j;
				first=false;
			}
			ind_end=j;

			temp=exp(temp);
			sum+=temp*tempvecNK_1[j];

		}
		IndJ_begin[i]=ind_begin;
		IndJ_end[i]=ind_end+1;

		if(sum>0)F[i]=1./sum;
		else F[i]=1;
	
	}
	vdLn(NK,F,tempvecNK_3);
	maxind=cblas_idamax(NK,tempvecNK_3,1);
	maxdif=fabs(tempvecNK_3[maxind]/k);

}


void do_sinkhorn_pure_axb(int k, double errortolerance)
{

	name+="Regparam"+to_string(k)+"_pure_axb";
	cout<<"Regularization parameter: "<<k<<endl<<"Expected Error tolerance:	"<<errortolerance<<endl;

	clock_t sink1;

	int i=0;
	int regvariable=multiplier*getk(NK_small);

	UseSmall(regvariable);

	while(regvariable<k)
	{
		sink1 = clock();

		Sinkhorn_axb(regvariable);
		absorbtion(regvariable);

		sink1= clock()-sink1;
		printf(" %d, regvariable %d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i,regvariable, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);

		regvariable+=pow(k,1./3.);

	}

	cout<<"reached regparam"<<endl;
	//one last iteration with final regularization parameter.

	i=0;
	while(maxdif>cap_treshold)
	{	
		sink1 = clock();
		Sinkhorn_axb(k);
		absorbtion(k);
		sink1= clock()-sink1;
		printf("%d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);
		if(i>cap_iteration)break;
	}

	cout<<"Number of iterations: "<<i<<endl<<"Last change in iterations: "<<maxdif<<endl;

		sink_to_Reflector_subtracted_axb(k);

}


void do_sinkhorn_subtracted_axb(int k, double errortolerance)
{

	name+="Regparam"+to_string(k)+"_subtracted_axb";
	cout<<"Regularization parameter: "<<k<<endl<<"Expected Error tolerance:	"<<errortolerance<<endl;

	clock_t sink1;

	int i=0;
	int regvariable=multiplier*getk(NK_small);

	UseSmall(regvariable);

	while(regvariable<k)
	{
		sink1 = clock();

		Sinkhorn_axb(regvariable);
		absorbtion(regvariable);

		sink1= clock()-sink1;
		printf(" %d, regvariable %d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i,regvariable, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);

		regvariable+=pow(k,1./3.);

	}

	cout<<"reached regparam"<<endl;
	//one last iteration with final regularization parameter.

	i=0;
	while(maxdif>cap_treshold)
	{	
		sink1 = clock();
		Sinkhorn_axb(k);
		absorbtion(k);
		sink1= clock()-sink1;
		printf("%d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);
		if(i>cap_iteration)break;
	}

	cout<<"Number of iterations: "<<i<<endl<<"Last change in iterations: "<<maxdif<<endl;



	cout<<endl<<"Started computing identity on first marginal"<<endl;

	i=0;
	regvariable=1;
	while(regvariable<k)
	{
		sink1 = clock();

		Sinkhorn_identity_F_axb(regvariable);
		absorbtion_f_id(regvariable);

		sink1= clock()-sink1;
		printf("%d,regvariable %d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i, regvariable, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);

		regvariable+=sqrt(k);

	}

	i=0;
	while(maxdif>cap_treshold)
	{
		sink1 = clock();
		Sinkhorn_identity_F_axb(k);
		absorbtion_f_id(k);
		sink1= clock()-sink1;
		printf("%d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);
		if(i>cap_iteration)break;

	}


	fill_n(IndId_begin,NK,0);
	fill_n(IndId_end,NK, NK);

	cout<<"Number of iterations: "<<i<<endl<<"Last change in iterations: "<<maxdif<<endl;

	cout<<endl<<"Started computing identity on second marginal"<<endl;

	i=0;
	regvariable=1;
	while(regvariable<k)
	{
		sink1 = clock();

		Sinkhorn_identity_G_axb(regvariable);
		absorbtion_g_id(regvariable);

		sink1= clock()-sink1;
		printf("%d,regvariable %d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i, regvariable, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);

		regvariable+=sqrt(k);
	}
	i=0;
	while(maxdif>cap_treshold)
	{
		sink1 = clock();
		Sinkhorn_identity_G_axb(k);
		absorbtion_g_id(k);
		sink1= clock()-sink1;
		printf("%d, time in clicks %ld, time in seconds %f with change %.*e \n", ++i, sink1, (float)sink1/CLOCKS_PER_SEC,DECIMAL_DIG,maxdif);
		if(i>cap_iteration)break;
	}
	cout<<"Number of iterations: "<<i<<endl<<"Last change in iterations: "<<maxdif<<endl;


	sink_to_Reflector_subtracted_axb(k);

}



void sink_to_Reflector_subtracted_axb(int k)
{
	Compute_TotalCost();

	//I do not understand how to normilize identities, as it seems they can have only one fixed value from received iterations, but they clearly need to be normalized


	double max =-1./0.;
	for(int i=0; i<NK; i++)
	{
		if(f_id[i]>max)max=f_id[i];
	}

	for(int i=0; i<NK; i++)
	{
		f_id[i]-=max;
		g_id[i]+=max;
	}

	max =-1./0.;
	for(int i=0; i<NK; i++)
	{
		if(f[i]>max)max=f[i];

	}

	for(int i=0; i<NK; i++)
	{
		f[i]-=max;
		g[i]+=max;
	}

	for(int i=0; i<NK; i++)
	{
		f[i]-=f_id[i];
		g[i]-=g_id[i];
	}


	vdExp(NK,f,R);
	
	for(int i=0; i<NK; i++)
	{
		for(int j=0; j<dim; j++)
		{
			Ref[i][j]=2*x[i][j]*R[i];
		}
	}



	//getting solution on regular grid


	clock_t transform=clock();


	for(int i=0; i<FinalGrid; i++)
	{
		double min=1./0.;
		for(int j=0; j<NK; j++)
		{
			if(g[j]>=1./0.)continue;
			double var=Cost_Func(x_regular[i],y[j])-g[j];
			if(var<min)min=var;
		}
		f_regular[i]=min;
	}

	for(int i=0; i<FinalGrid; i++)
	{
		for(int j=0; j<dim; j++)
		{
			Ref_regular[i][j]=2*x_regular[i][j]*exp(f_regular[i]);
		}
	}

	transform=clock()-transform;
	printf("Transform time in clicks %ld, time in seconds %f \n", transform, (float)transform/CLOCKS_PER_SEC);


	command=outputname+"/f_id_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(f_id);

	command=outputname+"/g_id_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(g_id);

	command=outputname+"/F_id_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(F_id);

	command=outputname+"/G_id_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(G_id);


	// command=outputname+"/f_regular_"+to_string(int(Reg_Param))+".txt";
	// freopen(command.c_str(),"w",stdout);
	// printf("%d  \n",FinalGrid);
	// print(f_regular,FinalGrid);

	// command=outputname+"/Ref_regular.txt";
	// freopen(command.c_str(),"w",stdout);
	// printf("%d  \n",FinalGrid);
	// print(Ref_regular,FinalGrid,dim);

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);

}

void absorbtion(int k)
{

	vdLn(NK,F,tempvecNK_1);	
	cblas_dscal(NK,1./(2*k),tempvecNK_1,1);
	vdAdd(NK,f,tempvecNK_1,tempvecNK_2);
	swap(f,tempvecNK_2);

	vdLn(NK,G,tempvecNK_1);
	cblas_dscal(NK,1./(2*k),tempvecNK_1,1);
	vdAdd(NK,g,tempvecNK_1,tempvecNK_2);
	swap(g,tempvecNK_2);

	fill_n(F,NK,1);
	fill_n(G,NK,1);

}

void absorbtion_f_id(int k)
{

	vdLn(NK,F_id,tempvecNK_1);
	cblas_dscal(NK,1./(2*k),tempvecNK_1,1);
	vdAdd(NK,f_id,tempvecNK_1,tempvecNK_2);
	swap(f_id,tempvecNK_2);

	fill_n(F_id,NK,1);

}

void absorbtion_g_id(int k)
{


	vdLn(NK,G_id,tempvecNK_1);
	cblas_dscal(NK,1./(2*k),tempvecNK_1,1);
	vdAdd(NK,g_id,tempvecNK_1,tempvecNK_2);
	swap(g_id,tempvecNK_2);

	fill_n(G_id,NK,1);
}


void GetCtransforms()
{
	Get_fc();
	Get_gc();

	for(int i=0; i<NK; i++)
	{
		for(int j=0; j<dim; j++)
		{
			Refc[i][j]=2*x[i][j]*exp(gc[i]);
		}
	}

	double dif_f[NK];
	double dif_g[NK];

	for(int i=0; i<NK; i++)
	{
		dif_f[i]=f[i]-gc[i];
		dif_g[i]=g[i]-fc[i];
	}

	command=outputname+"/Refc_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(Refc);

	command=outputname+"/fc_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(fc);

	command=outputname+"/gc_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(gc);

	command=outputname+"/dif_fc_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(dif_f);

	command=outputname+"/dif_gc_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(dif_g);

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);

	name+="_C";

}

void Get_fc()
{
	for(int j=0; j<NK; j++)
	{
		double min=1./0.;
		for(int i=0; i<NK; i++)
		{
			double var=Cost_Func(x[i],y[j])-f[i];
			if(var<min)min=var;
		}
		fc[j]=min;
	}
}

void Get_gc()
{
	for(int i=0; i<NK; i++)
	{
		double min=1./0.;
		for(int j=0; j<NK; j++)
		{
			if(g[j]>=1./0.)continue;
			double var=Cost_Func(x[i],y[j])-g[j];
			if(var<min)min=var;
		}
		gc[i]=min;
	}
}





void Printall_MY()
{

	command=outputname+"/x_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(x);


	command=outputname+"/y_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(y);

	command=outputname+"/Ref_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(Ref);

	command=outputname+"/R_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(R);

	command=outputname+"/p_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(p);

	command=outputname+"/q_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(q);

	command=outputname+"/F_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(F);

	command=outputname+"/G_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(G);

	command=outputname+"/f_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(f);

	command=outputname+"/g_MY.txt";
	freopen(command.c_str(),"w",stdout);
	printf("%d  \n",NK);
	print(g);

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);
}




void print(double vector [], int size)
{
	for(int i=0; i<size; i++)
	{
		printf("%.*e \n",DECIMAL_DIG,vector[i]);
	}
	printf("\n\n");
}

void print(int vector [], int size)
{
	for(int i=0; i<size; i++)
	{
		printf("%d \n",vector[i]);
	}
	printf("\n\n");
}


void print(double vector [][dim], int size1, int size2)
{
	for(int i=0; i<size1; i++)
	{
		for(int j=0; j<size2; j++)
		{
			printf("%.*e ",DECIMAL_DIG,vector[i][j]);

		}
		printf("\n");
	}
	printf("\n\n");
}


void UseSmall(int regparam)
{

	smallsinkhorn(regparam);

	for(int j=0; j<NK; j++)
	{
		double min=1./0.;
		for(int i=0; i<NK_small; i++)
		{
			if(f_small[i]>=1./0.)continue;
			double var=Cost_Func(x_small[i],y[j])-f_small[i];
			if(var<min)min=var;
		}
		g[j]=min;
	}

	
	for(int i=0; i<NK; i++)
	{
		double min=1./0.;
		for(int j=0; j<NK_small; j++)
		{
			if(g_small[j]>=1./0.)continue;
			double var=Cost_Func(x[i],y_small[j])-g_small[j];
			if(var<min)min=var;
		}
		f[i]=min;

	}

}

void Get_Original_Y(string outputname)
{
		command=outputname+"/Y_projected.txt";
	freopen(command.c_str(),"w",stdout);
	for(int i=0; i<NK; i++)
	{
		if(q[i]==0)continue;
		double X=y[i][0]/(1-y[i][2]);
		double Y=y[i][1]/(1-y[i][2]);
		
		double Z=Q(y[i]);

		printf("%.*e %.*e %.*e\n",DECIMAL_DIG,X,DECIMAL_DIG,Y,DECIMAL_DIG,Z);
	

	}


	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);
}



double Compute_TotalCost()
{

	double sum1 = 0;
	for(int i=0; i<NK; i++)
	{
		if(p[i]==0) continue;
		sum1+=p[i]*(f[i]-f_id[i]);
	}

	double sum2 = 0;
	for(int j=0; j<NK; j++)
	{
		if(q[j]==0) continue;
		sum2+=q[j]*(g[j]-g_id[j]);
	}

	cout<<endl<<sum1<<" + "<<sum2<<" = "<<sum1+sum2<<endl;
	return sum1+sum2;
}