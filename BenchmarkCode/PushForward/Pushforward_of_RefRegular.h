#include "../QuasiMonteCarlo/Generic_3D_logcost_MonteCarlo.h"

#include "../PushForward/PushForward_Cloud_128.h"

/*


For all x-montecarlo
	find patch where x belongs
	compute bi-linear approximation
	compute direction of the reflection & density
	save projection of the direction and density



*/

double coef[4];
vector <double> Do_RegularPush(double x,double y,double z);

int binsearch(int a, int b, double val, double vect[]);
void IntCoef(int i, int j);
double rightrootofquadraticequation(double a, double b, double c, double average);



void Pushforward_Ref_regular(string outputname)
{
	string command=outputname+"/Y_Pushed_projected.txt";
	freopen(command.c_str(),"w",stdout);

	for (int ind=0; ind<Push_Cloud_Size; ind++)
	{
		vector <double> vec=Do_RegularPush(Push_Cloud[ind][0],Push_Cloud[ind][1],Push_Cloud[ind][2]);
		if(vec[3]==0) continue;
		printf("%.*e %.*e %.*e \n",	DECIMAL_DIG,vec[0]/(1-vec[2]),
									DECIMAL_DIG,vec[1]/(1-vec[2]),
									DECIMAL_DIG,vec[3]);		
	}	

	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);


	for (int ind=0; ind<NK; ind++)
	{
		vector <double> vec=Do_RegularPush(x[ind][0],x[ind][1],x[ind][2]);
		Y_pushed_projected[ind][0]=vec[0]/(1-vec[2]);
		Y_pushed_projected[ind][1]=vec[1]/(1-vec[2]);
		Y_pushed_projected[ind][2]=vec[3];

	}

	command=outputname+"/Y_Pushed_projected_owndisc.txt";
	freopen(command.c_str(),"w",stdout);


	for (int ind=0; ind<NK; ind++)
	{
		if(Y_pushed_projected[ind][2]==0) continue;
		printf("%.*e %.*e %.*e \n",	DECIMAL_DIG,Y_pushed_projected[ind][0],
									DECIMAL_DIG,Y_pushed_projected[ind][1],
									DECIMAL_DIG,Y_pushed_projected[ind][2]);		
	}


	command=outputname+"/X_projected.txt";
	freopen(command.c_str(),"w",stdout);

;
	for (int ind=0; ind<NK; ind++)
	{
		if(P(x[ind])==0) continue;
		double X=x[ind][0]/(1+x[ind][2]);
		double Y=x[ind][1]/(1+x[ind][2]);
		printf("%.*e %.*e %.*e \n",	DECIMAL_DIG,X,
									DECIMAL_DIG,Y,
									DECIMAL_DIG,P(x[ind]));		
	}


	command=outputname+"/log.txt";
	freopen(command.c_str(),"a",stdout);
}

vector <double> Do_RegularPush(double x,double y,double z)
{
	double X=x/(1+z);
	double Y=y/(1+z);
	// find corresponding patch to be done

	int i=binsearch(0,FinalGridResolution, X, Regular_side);
	int j=binsearch(0,FinalGridResolution, Y, Regular_side);


	IntCoef(i,j);  // generate interpolating plane (indexes  should come from previous step)

	//computing scale for vector by solving equation coming from reflector interpolation
	//hz=a0+a1hx+a2hy+a3hxhy

	double averagelength= ( L2_norm(Ref_regular[i*FinalGridResolution+j])+L2_norm(Ref_regular[i*FinalGridResolution+j+1])
							+ L2_norm(Ref_regular[(i+1)*FinalGridResolution+j])+L2_norm(Ref_regular[(i+1)*FinalGridResolution+j+1]) ) /4;

	double h=rightrootofquadraticequation(coef[3]*x*y, coef[1]*x+coef[2]*y-z, coef[0], averagelength);

	double tempvar1=(coef[1]+coef[3]*h*y);//-n1
	double tempvar2=(coef[2]+coef[3]*h*x);//-n2
	double tempvar3=2*(coef[0]/h-coef[3]*h*x*y)
					/ (tempvar1*tempvar1+tempvar2*tempvar2+1);
					//(2*<a,n>/normalsnorm)

	vector <double> vec (4);

	vec[0]=x + tempvar1*tempvar3; // x1-n1*2<x,n>
	vec[1]=y + tempvar2*tempvar3; // x2-n2*2<x,n>
	vec[2]=z - tempvar3; //x3-n3*2<x,n> but n3=1



	double * Temp=new double [3];
	Temp[0]=x;
	Temp[1]=y;
	Temp[2]=z;

	vec[3]=P(Temp);//*4/((1+X*X+Y*Y)*(1+X*X+Y*Y));  //density from the sphere times the jacobian of stereographic projection.


	return vec;
	//Return values of reflected direction vector (a,b,c)
	//and value of density of incoming ray
}
































int binsearch(int a, int b, double val, double vect[])
{
	int mid=(a+b)/2;
	if(mid==a) return a;  //if interval already found, return interval index

	if(val<vect[mid]) return binsearch(a, mid, val, vect);
	else return binsearch (mid, b, val, vect);
}



void IntCoef(int i, int j)
{
	int k=FinalGridResolution;
	/*

	Bi-linear interpolation on 4 points in space: coef[0]+coef[1]x+coef[2]y+cpef[3]xy=f(x,y)

	equations are: 
	coef[0]+coef[1]x1+coef[2]y1+cpef[3]x1y1=z1;
	coef[0]+coef[1]x2+coef[2]y2+cpef[3]x2y2=z2;
	coef[0]+coef[1]x3+coef[2]y3+cpef[3]x3y3=z3;
	coef[0]+coef[1]x4+coef[2]y4+cpef[3]x4y4=z4;


	*/

	//gauss elimination written by hand (to do it faster hopefully....)

	//indexation taken as 1=ij 2=i+1j 3=ij+1 4=i+1j+1

	double xdif2=Ref_regular[(i+1)*k+j][0]-Ref_regular[i*k+j][0];
	double xdif3=Ref_regular[(i)*k+j+1][0]-Ref_regular[i*k+j][0];
	double xdif4=Ref_regular[(i+1)*k+j+1][0]-Ref_regular[i*k+j][0];

	double ydif2=Ref_regular[(i+1)*k+j][1]-Ref_regular[i*k+j][1];
	double ydif3=Ref_regular[(i)*k+j+1][1]-Ref_regular[i*k+j][1];
	double ydif4=Ref_regular[(i+1)*k+j+1][1]-Ref_regular[i*k+j][1];

	double fdif2=Ref_regular[(i+1)*k+j][2]-Ref_regular[i*k+j][2];
	double fdif3=Ref_regular[(i)*k+j+1][2]-Ref_regular[i*k+j][2];
	double fdif4=Ref_regular[(i+1)*k+j+1][2]-Ref_regular[i*k+j][2];

	double pdif2=Ref_regular[(i+1)*k+j][0]*Ref_regular[(i+1)*k+j][1]-Ref_regular[i*k+j][0]*Ref_regular[i*k+j][1];
	double pdif3=Ref_regular[(i)*k+j+1][0]*Ref_regular[(i)*k+j+1][1]-Ref_regular[i*k+j][0]*Ref_regular[i*k+j][1];
	double pdif4=Ref_regular[(i+1)*k+j+1][0]*Ref_regular[(i+1)*k+j+1][1]-Ref_regular[i*k+j][0]*Ref_regular[i*k+j][1];

	double by3=ydif3-ydif2*xdif3/xdif2;
	double by4=ydif4-ydif2*xdif4/xdif2;

	double bp3=pdif3-pdif2*xdif3/xdif2;
	double bp4=pdif4-pdif2*xdif4/xdif2;	

	double bf3=fdif3-fdif2*xdif3/xdif2;
	double bf4=fdif4-fdif2*xdif4/xdif2;	

	coef[3]=(bf4-bf3*by4/by3)/(bp4-bp3*by4/by3);
	coef[2]=(bf3-coef[3]*bp3)/by3;
	coef[1]=(fdif2-coef[3]*pdif2-coef[2]*ydif2)/xdif2;
	coef[0]=Ref_regular[i*k+j][2]-coef[3]*Ref_regular[i*k+j][0]*Ref_regular[i*k+j][1]-coef[2]*Ref_regular[i*k+j][1]-coef[1]*Ref_regular[i*k+j][0];


}




double rightrootofquadraticequation(double a, double b, double c, double average)
{
	double root1,root2;
	//accroding to some one page paper from MIT, computing roots of quadratic equation might have rounding error when b is bigger then a*c
	//even thou it's not the case here, there isn't much effort in computing it by their suggested way, so why not

	double diskr=b*b-4*a*c;
	if(diskr<0) 
	{
		cout<<"Error!, quadratic equation doesn't have real solution! It shouldn't have happened according to model!"<<endl;
		return -1;
	}
	diskr=sqrt(diskr);

	if(b<0)
	{
		root1=(2*c)/(-b+diskr);
		root2=(-b+diskr)/(2*a);
	}
	else 
	{
		root1= (-b-diskr)/(2*a);
		root2= (2*c)/(-b-diskr);
	}

	if(fabs(root1-average)<fabs(root2-average))	return root1;
	else return root2;
	
}