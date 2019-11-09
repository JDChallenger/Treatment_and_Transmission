
#include "parameters.h"
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

//////////////////////////////////////////////////
/* Define RK4 functions */
//////////////////////////////////////////////////
double f1AM(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = -(theta->kaAM) * x;
return temp;
}

double f2AM(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->kaAM) * x -((theta->k23AM) + ((theta->CLAM) / (theta->VCAM)))* y;
return temp;
}

double f3AM(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->k23AM) * y - ((theta->CLmetAM) / (theta->VMAM))*z;
return temp;
}

double f1L(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = -(theta->kaL) * x;
return temp;
}

double f2L(double x,double y, double z , params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->kaL) * x -((theta->k23L) + ((theta->CLL) / (theta->VCL)))* y;
return temp;
}

double f3L(double x,double y, double z, params* theta)//(x,y,z)=(gut,central,metabolite)
{
double temp;
temp = (theta->k23L) * y - ((theta->CLmetL) / (theta->VML))*z;
return temp;
}


//And now for Gametocytes. Simple structure: in and out
//Used for untreated infections
double gf(double x, double y, double nux, double nuy) //In / out / in / out . Could put nu in a structure
{
double temp;
temp = nux * x - nuy * y;
return temp;
}

//Used for treated infections
double gfPD(double x, double y, double nux, double nuy, int a, double AMconc, double DHAconc, double LMFconc, paramsPD2* thetaPD2) 
{
double temp;
temp = nux * x - nuy * y - ((thetaPD2->Gkmax_AM[a] * AMconc /(AMconc + thetaPD2->c50_AM)) + 
	(thetaPD2->Gkmax_AM[a] * DHAconc /(DHAconc + thetaPD2->c50_AM)) + (thetaPD2->Gkmax_L[a] * LMFconc /(LMFconc + thetaPD2->c50_L)) )*y;
return temp;
}

//Bradley model, using (inferred) male & female gametocyte populations 
int infectivityB2(double xg)
{
	//double temp, xgM, xgF;//Don't need denominator for A3=0
	double G10 = log10(xg);
	int guess = (4 + G10)*200;
	//temp = (1-exp(-thetaIB->tau * xgM))* thetaIB->A1 * pow(xgF,thetaIB->A2);
	if(guess>1562)
		guess = 1562;
	if(xg<0.00011)
		guess = 1;
	return guess;
}

//Implementation of the analytical solution for the Artemether PK equation
vector<double> ART(double Time, double Dose, params* theta){
	vector<double> temp(3);
	
	double c1 = - (Dose * theta->k23AM * theta->kaAM * theta->VCAM * theta->VMAM * theta->VMAM)/
	((-theta->CLmetAM + theta->kaAM*theta->VMAM)*
		(-theta->CLmetAM * theta->VCAM + theta->CLAM * theta->VMAM + theta->k23AM * theta->VCAM * theta->VMAM));
	double c2 = (Dose * theta->k23AM * theta->kaAM * theta->VCAM * theta->VMAM)/
	((-theta->CLmetAM + theta->kaAM*theta->VMAM)*
		(theta->CLAM + theta->k23AM * theta->VCAM - theta->kaAM * theta->VCAM));
	double c3 = -(Dose * theta->k23AM * theta->kaAM * theta->VCAM * theta->VCAM * theta->VMAM)/
	((-theta->CLmetAM * theta->VCAM + theta->CLAM * theta->VMAM + theta->k23AM * theta->VCAM * theta->VMAM)*
		(theta->CLAM + theta->k23AM * theta->VCAM - theta->kaAM * theta->VCAM));

	double eig1 = -theta->CLmetAM / theta->VMAM;
	double eig2 = -theta->kaAM;
	double eig3 = -(theta->CLAM + theta->k23AM * theta->VCAM)/(theta->VCAM);

	double eigV1[3] = {0,0,1};
	double eigV2[3] = 
		{((theta->CLAM + theta->k23AM * theta->VCAM - theta->kaAM * theta->VCAM)*
			(-theta->CLmetAM + theta->kaAM * theta->VMAM))/(theta->k23AM * theta->kaAM * theta->VCAM * theta->VMAM),
		-(theta->CLmetAM - theta->kaAM * theta->VMAM)/(theta->k23AM * theta->VMAM),1};
	double eigV3[3] = {0,
		-(theta->CLmetAM * theta->VCAM - theta->CLAM * theta->VMAM - theta->k23AM * theta->VCAM * theta->VMAM)/
		(theta->k23AM * theta->VCAM * theta->VMAM),1};

	//Eigenvector #3 always has a factor of -1 which we need to remove
	temp[0] = c1 * exp(eig1*Time) * eigV1[0] + c2 * exp(eig2*Time) * eigV2[0] + c3 * exp(eig3*Time) * eigV3[0];
	temp[1] = c1 * exp(eig1*Time) * eigV1[1] + c2 * exp(eig2*Time) * eigV2[1] + c3 * exp(eig3*Time) * eigV3[1];
	temp[2] = abs(c1 * exp(eig1*Time) * eigV1[2] + c2 * exp(eig2*Time) * eigV2[2] + c3 * exp(eig3*Time) * eigV3[2]); 

	return temp;
}


