#include <iostream>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <random>
#include <cstdlib>
#include "parameters.h"
//#include <ctime>
#include <cstring>
#include <vector>
#include <chrono>

using namespace std;

int main () {

//Name of output file, in which to store the simulation
const char output_File[] = "WHM_Gametocytes_Results.txt";

//****** 1. Read in data file which determines the infectivity model ****** 

//Input library file For the infectivity model due to Bradley and co-workers
const char infilename[] = "Bradley_Infectivity.dat";
int mx = 1563;
int max = mx*5;

double vv;
vector<double>vec;

ifstream in(infilename);
	if(!in){
	cerr <<"Failed to open input file "<< infilename << endl;
exit(1);
	}
	int ii=0;
while(in){
  if(ii >= max) break;
  	in >>vv;
  	vec.push_back(vv);
  	ii++;
	}
	in.clear();
	in.close();

//cout<<"Vector vec contains "<< vec.size() <<" entries. "<<vec[1]<<endl;
double lib_gtot[mx];
double lib_gM[mx];
double lib_gF[mx];
double lib_ratio[mx]; //Note: this is M/F
double lib_inf[mx]; //Note: this is M/F
for(int i=0;i<mx;i++){
	lib_gtot[i] = vec[5*i];
	lib_gF[i] = vec[5*i + 1];
	lib_gM[i] = vec[5*i + 2];
	lib_ratio[i] = vec[5*i + 3];
	lib_inf[i] = vec[5*i + 4];
}

//******  2. Define all parameter values for this simulation ****** 

//Generate enough correlated random numbers for one run here?
vector<double> R1(450,0.0); //These will be uncorrelated random numbers
vector<double> R ; //Random numbers correlated in time

//construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

normal_distribution<double> distribution (0,1);

uniform_real_distribution<double> distribution2(0.0,1.0);
//cout <<"Uniform random variable: "<<distribution2(generator) <<endl;

for(int j=0;j<450;j++){
	double number = distribution(generator);
	R1[j] = number;
}

//Using uncorrelating random numbers above, generate random numbers correlated in time
double sum=0;
double bray=0;
for(int j=0;j<450;j++){
	sum=0;
	for(int j1=1;j1<j+1;j1++){
		sum+=R1[j1]*pow(f,j-j1);
	}
	bray = pow(f,j)*R1[0]+sqrt(1-pow(f,2))*sum;
	if(mu + sigma * bray>1 && mu + sigma * bray<35){
  		R.push_back (mu + sigma * bray);
	}
}
//cout<<"Check length of R: "<<R.size() <<endl;

//Now generate the random numbers needed for the innate & general-adaptive immunities
double Pms;
double Pcs;

double meanLN = 4.79*log(10);
double sigmaLN =1.2;

double rand1 = meanLN + sigmaLN * distribution(generator);
double rand2 = exp(rand1);

while(rand2 > pow(10,5.5)){
	rand1 = meanLN + sigmaLN * distribution(generator);
	rand2 = exp(rand1);
}

Pcs = kc * rand2;
//cout<<"Pcs: " << Pcs <<endl;

Pms = km * log(1-(1/gompcst2)*(log(1-distribution2(generator))))*(1/gompcst1);
//cout<< "Pms: " << Pms <<endl;

/* Declare vector to store parasitaemia, and give initial condition at t=0 */
vector<double> P(400,0.0);
P[0]=0.1;//initial condition

vector<double> PC2(400,0.0); //Auxilliary vector, for calculating adaptive immune response
double PC;
double PCsum;
double Pv=0;

double Sc = 1;
double Sm = 1;
double Sv = 1;
int upper, lower;

//Choose a value of dt (often called h in RK4)
double dt = 0.25;
int delta_t=1/dt;
int TL = (1/dt) * 24 * 800;// For 400 generations' worth
cout <<"For 800 days we require "<< TL <<" increments" <<endl;

//Also useful to have the number of intervals for one & two days (one generation)
int oneday = (1/dt) * 24;
int twoday = (1/dt) * 48;

int L = 3; //Number of aux. compartments at each gametocyte life stage.

vector<double> aux(TL,0.0);
vector<vector<double> > G(5 * L,aux);
cout<<"Dimensions of G: " << G.size() <<endl;
aux.clear();

//Initialise G1
G[0][0] = 0.0;

//Gametocyte parameters: Draw values for alpha, delta and delta5 from relevant distributions.
double alpha = exp( -alpha_mean + alpha_sd * distribution(generator) );
while(alpha>0.3){
	alpha = exp( -alpha_mean + alpha_sd * distribution(generator) );
}
cout<<"Sexual commitment rate: "<< alpha <<endl;

double delta = delta_mean + delta_sd * distribution(generator) ;
while(delta<24||delta>72){
	delta = delta_mean + delta_sd * distribution(generator) ;
}
cout<<"Delta: "<< delta <<endl;

double delta5 = delta5_mean + delta5_sd * distribution(generator);//exp( 4.98 + 0.765 * distribution(generator) );
while(delta5<33.6||delta5>384){
	delta5 = delta5_mean + delta5_sd * distribution(generator);//exp( -4.15 + 2.093 * distribution(generator) );
}
cout<<"Delta5: "<< delta5 <<endl;

double nu = (1/delta) * L; //rate of transitions through gametocyte compartments
double d5 = (1/delta5) * L; //death rate of mature gametocytes

double nuv[5] = {nu,nu,nu,nu,d5};
double fG1[L*5] , fG2[L*5] , fG3[L*5], fG4[L*5]; //arrays for Runge-Kutta routine

int uu1;
int uu2;

int AS = 0;
int CF = 0;
int fevk = 0;

//****** 	 3. Run the model	****** 

for(int k = 0; k < TL - 1; k++){

	fG1[0] = dt * gf(0,G[0][k],0,nuv[0]);
	fG2[0] = dt * gf(0,G[0][k] + 0.5 * fG1[0],0,nuv[0]);
	fG3[0] = dt * gf(0,G[0][k] + 0.5 * fG2[0],0,nuv[0]);
	fG4[0] = dt * gf(0,G[0][k] + fG3[0],0,nuv[0]);
	//RK4 step. Combine these four fns together
	G[0][k+1] = G[0][k] + (fG1[0] + 2 * fG2[0] + 2 * fG3[0] +fG4[0])/6;
	for(int a=1; a < L * 5; a++){
		uu1 = (a-1)/L;
		uu2 = a/L;
		fG1[a] = dt * (gf(G[a-1][k],G[a][k],nuv[uu1],nuv[uu2]));
		fG2[a] = dt * (gf(G[a-1][k] + 0.5 * fG1[a-1],G[a][k] + 0.5 * fG1[a],nuv[uu1],nuv[uu2]));
		fG3[a] = dt * (gf(G[a-1][k] + 0.5 * fG2[a-1],G[a][k] + 0.5 * fG2[a],nuv[uu1],nuv[uu2]));
		fG4[a] = dt * (gf(G[a-1][k] + fG3[a-1],G[a][k] + fG3[a],nuv[uu1],nuv[uu2]));
		//RK4 step. Combine these four fns together
		G[a][k+1] = G[a][k] + (fG1[a] + 2 * fG2[a] + 2 * fG3[a] +fG4[a])/6;
	}

	//Every 48 hours, update the asexual parasitaemia model (see 2017 Nat Comms paper)
	if(k%(48 * delta_t)==twoday-1 && AS==0){ 
		int k1 = ((k+1) * dt /48 )-1;

		//Innate immune response
		Sc = 1/(1+pow((P[k1]/Pcs),kaC));

		//Effective var-specific adaptive immune response
		//First determine the limits on the sum
		if(k1>3){
			lower = round((k1+1-4)*pow(lambda , k1-4+1))-1;
			upper = k1 - 4 +1;
			Pv=0;
			for(int k2=lower;k2<upper;k2++){
				Pv += P[k2];
			}
			Sv=1/(1+pow((Pv/Pvs),kaV));
		}
		// General-adaptive immune response
		if(P[k1]>C){
			PC=C;
		}
		else{
			PC=P[k1];
		}

		PC2[k1]=PC;
		
		if(k1>3){
			PCsum=0;
			for(int k2=0;k2<k1-4+1;k2++){
				PCsum+=PC2[k2];
			}
			Sm = beta + (1-beta)/(1+pow((PCsum/Pms),kaM));
		}

		P[k1+1] =  R[k1] * Sc * Sm * Sv * P[k1];
		G[0][k+1] += alpha * R[k1] * Sc * Sm * Sv * P[k1];
		//cout<<"On Day "<< 2*k1<<", " << P[k1]<<" "<<Sc<<" "<<Sm<<" "<<Sv<<" "<< G[0][k+1]<<endl;

		//Has the asexual parasitaemia been cleared?
		if(P[k1+1]<pow(10,-5)){
			AS=1;
		}
	}

	//Check if there are parasites left- if not, break out of the loop
	if(AS==1){
		double gam_mature = 0.0;
		for(int ii=0;ii<L;ii++){
			gam_mature += G[4*L + ii][k];
		}
		if(AS>0 && gam_mature < pow(10,-5)){
			cout<< "Last mature gametocytes have gone" <<endl;
			break;
		}
	}

}//end of 'dt' loop

//******  4. Simulation finished, prepare data for saving to file. ****** 

vector<double> GT1(TL,0.0);
vector<double> GT2(TL,0.0);
vector<double> GT3(TL,0.0);
vector<double> GT4(TL,0.0);
vector<double> GT5(TL,0.0);

for(int k=0;k<TL;k++){
	for(int i=0;i<L;i++){
		GT1[k] += G[i][k];
		GT2[k] += G[L + i][k];
		GT3[k] += G[2*L + i][k];
		GT4[k] += G[3*L + i][k];
		GT5[k] += G[4*L + i][k];
	}
}
G.clear();

//Convert gametocyte density to infectivity (using model due to Bradley et al.)
vector<double> P_infect(TL,0.0);
for(int k=0;k<TL;k++){
	P_infect[k] = lib_inf[infectivityB2(GT5[k])];
}
//Now calculate AUIC
double auic = 0.0;
for(int k=0;k<TL-1;k++){
	auic += (dt/24)*0.5*(P_infect[k]+P_infect[k+1]);
}
cout<<"Area under the infectivity curve: "<< auic <<endl;

cout<<"Now writing results to output file"<<endl; 
 	{
	  ofstream out(output_File);
	  if(!out){
	    cerr <<"Failed to open output file"<< output_File << endl;
	    exit(1);
	  }
	  for(int f3=0;f3<400;f3++){
	    out << 2*f3 <<"\t"<<P[f3]<<"\t"<<GT5[f3*twoday]<<"\t"<< P_infect[f3*twoday] <<endl;
	  }
	    out.close();
	}

} //end of main
