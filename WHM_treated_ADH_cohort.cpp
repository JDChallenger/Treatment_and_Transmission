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

//Input library file for Bradley model
const char infilename1[] = "Bradley_Infectivity.dat";
//Can this stuff go elsewhere??
int mx = 1563;
int max = mx*5;//How many elements in file?

double vv1;
vector<double>vec1;

ifstream in1(infilename1);
if(!in1){
	cerr <<"Failed to open input file "<< infilename1 << endl;
	exit(1);
}
int ii=0;
while(in1){
  if(ii >= max) break;
  	in1 >>vv1;
	vec1.push_back(vv1);
	ii++;
}
in1.clear();
in1.close();

cout<<"Vector vec contains "<< vec1.size() <<" entries. "<<vec1[1]<<endl;
double lib_gtot[mx];
double lib_gM[mx];
double lib_gF[mx];
double lib_ratio[mx];//Note: this is M/F
double lib_inf[mx];//Note: this is M/F
for(int i=0;i<mx;i++){
	lib_gtot[i] = vec1[5*i];
	lib_gF[i] = vec1[5*i + 1];
	lib_gM[i] = vec1[5*i + 2];
	lib_ratio[i] = vec1[5*i + 3];
	lib_inf[i] = vec1[5*i + 4];
}

//Output file for results
const char output_File1[] = "WHM_Results_para_PK_ADH.txt";
//const char output_File1[] = "WHM_Results_para_gam_baseline.txt";
const char output_File2[] = "WHM_Results_gam_PK_ADH.txt";
const char output_File3[] = "WHM_Results_AM_LMF_ADH.txt";

ofstream out1(output_File1);
ofstream out2(output_File2);
ofstream out3(output_File3);

//Choose a value of dt (often called h in RK4)
double dt = 0.2;
int delta_t=1/dt;
int TL = (1/dt) * 24 * 400;// For 400 generations' worth
cout <<"For 400 days we require "<< TL <<" increments" <<endl;

//Also useful to have the number of intervals for one & two days (one generation)
int oneday = (1/dt) * 24;
int twoday = (1/dt) * 48;

//How many patients?
int N=659;
//How many times to simulate each adherence profile?
int NN = 10;
//If infection recrudesces, with what probability is the patient retreated?
double prob_retreat = 0.333;

/* START WITH READING IN ADHERENCE DATA */

//Input file for adherence	
const char infilename2[] = "myadh.txt";

/* Read in the adherence data */
int dosing[659];//variable for pills per dose. One int for each patient
double timings[659][24];//timings for each pill. Up to 24, file is 'padded' with '-1's

int max2=(N*25);//How many elements in file?

double vv;
vector<double>vec;

ifstream in2(infilename2);
  if(!in2){
    cerr <<"Failed to open input file "<< infilename2 << endl;
    exit(1);
  }
  int ii2 = 0;
    while(in2){
      if(ii2 >= max2) break;
      in2 >>vv;
      vec.push_back(vv);
      ii2++;
    }
  in2.clear();
  in2.close();
  
cout<<"Vector vec contains "<< vec.size() <<" entries. "<<vec[1]<<endl;

for(int i=0; i<N; i++){

	dosing[i]=vec[i*25];

	for(int j=0;j<24;j++){
		timings[i][j] = vec[(25*i) + j + 1];
	}

}
vec.clear();

vector<double>AUICvec(659,0.0);
for(int GG=0;GG<N;GG++){
//Instead of using the recommended dose timings, select an adherence profile to use
//int GG=G;// Make this clearer for user

cout<<"Patient # "<<GG;
for(int h=0;h<18;h++){
	cout<<" "<<timings[GG][h]<<" ";
}
cout<<" "<<endl;

int d24[24]={0};

//When did this patient take their pills?
int f_last = 0;//Record final dose time (in increments)
for(int k=0;k<24;k++){
	d24[k] = (1/dt) * timings[GG][k];
	if(d24[k]>f_last){
		f_last = d24[k];
	}
}
/*for(int h=0;h<18;h++){
	cout<<" "<<d24[h]<<" ";
}
cout<<" "<<endl;*/

vector<int> d11;
vector<int> dm11;
d11.push_back(0);
dm11.push_back(1);
//cout<<d11[0]<<endl;
int nd=0;
for(int h=1;h<24;h++){
	if(d24[h]<0)
		break;
	if(d24[h]==d24[h-1]){
		dm11[nd]+=1;
	}
	else{
		nd++;
		d11.push_back(d24[h]);
		dm11.push_back(1);
	}
	
}
int l1=d11.size();
int l2=dm11.size();

/*for(int h=0;h<l1;h++){
	cout<<d11[h]<<" ";
}
cout<<" "<<endl;
for(int h=0;h<l2;h++){
	cout<<dm11[h]<<" ";
}
cout<<" "<<endl;*/

/*Ok separate to the considerations needed for the ODE soln, we must decide if a dosing schedule
should be discarded due to overdosing.	*/

int TABLETS = vec[25*GG]; //No. of pills that should be in a dose
int dosenum = 0; //Start from one or zero?
vector<double> dosetimes;
dosetimes.push_back(0.0);
vector<int> doseSize;
doseSize.push_back(1);

for(int h=0;h<23;h++){
	if(timings[GG][h+1]<0)
		break;
	if(timings[GG][h+1] <= 0.5 + timings[GG][h]){
		doseSize[dosenum]+=1;
	}
	else{
		dosenum++;
		doseSize.push_back(1);
		dosetimes.push_back(timings[GG][h+1]);
	}
}
cout<<"START OF OVERDOSING ASSESSMENT" <<endl;
cout<<" " <<endl;
cout<<"No. of doses: "<< dosenum+1 <<endl;
cout<<"Timings of doses:";
for(int h=0;h<dosenum+1;h++){
	cout<<" "<<dosetimes[h];
}
cout<<" " <<endl;
cout<<"No. of pills in each dose:";
for(int h=0;h<dosenum+1;h++){
	cout<<" "<<doseSize[h];
}
cout<<" " <<endl;

int c1=0;
for(int h=0;h<dosenum+1;h++){
	if(doseSize[h]>TABLETS){
		cout<<"PANIC" <<endl;
		c1++;
		break;
	}
}
if(c1==0)
	cout<<"DON'T PANIC"<<endl;
cout<<" " <<endl;
cout<<"END OF OVERDOSING ASSESSMENT" <<endl;


//construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator (seed);

normal_distribution<double> distribution (0,1);

uniform_real_distribution<double> distribution2(0.0,1.0);
//cout <<"Uniform random variable: "<<distribution2(generator) <<endl;

uniform_real_distribution<double> distribution3(-3.69,-0.0);//Upper end changed back

int L = 3; //Number of aux. compartments at each gametocyte life stage.

for(int q=0;q<NN;q++){
	//Generate enough correlated random numbers for one run here? Inside q loop or no?
	vector<double> R1(250,0.0); //These will be uncorrelated random numbers
	vector<double> R ; //Random numbers correlated in time

	for(int j=0;j<250;j++){
		double number = distribution(generator);
		R1[j] = number;
	}

	//cout<<"Test "<<pow(f,0)<<endl;
	double sum=0;
	double bray=0;
	for(int j=0;j<250;j++){
		sum=0;
		for(int j1=1;j1<j+1;j1++){
			sum+=R1[j1]*pow(f,j-j1);
		}
	bray = pow(f,j)*R1[0]+sqrt(1-pow(f,2))*sum;
	//cout<< mu + sigma * bray <<endl;
	if(mu + sigma * bray>1 && mu + sigma * bray<35){
  		R.push_back (mu + sigma * bray);
	}
	}
	//cout<<"Check length of R: "<<R.size() <<endl;

	//Now generate the random numbers needed for the innate & general-adaptive immunities
	double Pms;
	double Pcs;

	double meanLN = 4.79*log(10);
	double sigmaLN = 1.2;

	double rand1 = meanLN + sigmaLN * distribution(generator);
	double rand2 = exp(rand1);
	//cout<< rand2 <<" upper truncation: "<<pow(10,5.5)<<endl;

	while(rand2 > pow(10,5.5)){
		rand1 = meanLN + sigmaLN * distribution(generator);
		rand2 = exp(rand1);
		//cout<<"Pcs Replaced"<<endl;
	}

	Pcs = kc * rand2;
	//cout<<"Pcs: " << Pcs <<endl;

	Pms = km * log(1-(1/gompcst2)*(log(1-distribution2(generator))))*(1/gompcst1);
	//cout<< "Pms: " << Pms <<endl;

	vector<double> PC2(400,0.0);
	double PC;
	double PCsum;
	double Pv=0;

	double Sc = 1;
	double Sm = 1;
	double Sv = 1;
	int upper, lower;


	vector<double> aux(TL,0.0);
	//cout<<"Dimensions of aux: " << aux.size() <<endl;
	vector<vector<double> > G(5 * L,aux);
	//cout<<"Dimensions of G: " << G.size() <<endl;

	aux.clear();

	//Initialise G1
	G[0][0] = 0.0;

	//Draw values for alpha, delta and delta5 from relevant distributions. Put these in a struct?
	//Sexual commitment rate
	double alpha = exp( -alpha_mean + alpha_sd * distribution(generator) );
	while(alpha>0.3){
		alpha = exp( -alpha_mean + alpha_sd * distribution(generator) );
	}
	//cout<<"Sexual commitment rate: "<< alpha <<endl;

	double delta = delta_mean + delta_sd * distribution(generator) ;
	while(delta<24||delta>72){
		delta = delta_mean + delta_sd * distribution(generator) ;
	}
	//cout<<"Delta: "<< delta <<endl;

	double delta5 = delta5_mean + delta5_sd * distribution(generator);//exp( 4.98 + 0.765 * distribution(generator) );
	while(delta5<33.6||delta5>384){
		delta5 = delta5_mean + delta5_sd * distribution(generator);//exp( -4.15 + 2.093 * distribution(generator) );
	}
	//cout<<"Delta5: "<< delta5 <<endl;

	double nu = (1/delta) * L; //rate of transitions through gametocyte compartments
	//double nu1 = (1/(delta +24)) * L; //Stage I gametocytes of longer duration? Drug action on IIa and IIb the same?
	double d5 = (1/delta5) * L; //death rate of mature gametocytes

	double nuv[5] = {nu,nu,nu,nu,d5};
	double fG1[L*5] , fG2[L*5] , fG3[L*5], fG4[L*5];

	//Default is no drugs (need fever as the trigger for seeking treatment), so add DELTA
	int DELTA = TL + 5;
	int DELTA_old = 0;
	//cout<<"Check timings: "<<D1 <<" "<<D2 <<" "<<D3 <<" "<<D4 <<" "<<D5 <<" "<<D6 <<endl;

	double DRUGL;
	//double DRUGAM;

	double GUTAM[TL];//RESET THE WHOLE OF GUT EACH RUN
	double CENTRALAM[TL];
	double METABOLITEAM[TL];
	GUTAM[0] = 0.0;
	CENTRALAM[0] = 0.0;
	METABOLITEAM[0] = 0.0;

	double GUTL[TL];//RESET THE WHOLE OF GUT EACH RUN
	double CENTRALL[TL];
	double METABOLITEL[TL];
	//GUTL[0] = 0.0;
	CENTRALL[0] = 0.0;
	METABOLITEL[0] = 0.0;

	//Then, use these to apply the units, and do the killing
	double CENTRALAMx[TL];
	double METABOLITEAMx[TL];
	CENTRALAMx[0] = 0.0;
	METABOLITEAMx[0] = 0.0;

	double CENTRALLx[TL];
	double METABOLITELx[TL];
	CENTRALLx[0] = 0.0;
	METABOLITELx[0] = 0.0;

	double r1, r2, r3, r4, r5, r6;
	//Now scale each one with the appropriate standard deviation
	r1 = 0.44 * distribution(generator); // CLAM
	r2 = 1.19 * distribution(generator); // kaAM
	r3 = 0.68 * distribution(generator); // k23AM
	r4 = 0.38 * distribution(generator); // CLL
	r5 = 0.6 * distribution(generator); // VCL
	r6 = 0.38 * distribution(generator); // k23L

	//////////////////////////////////////////////////
	/* Generate PK parameters in the struct. Put variables in the header? */
	//////////////////////////////////////////////////
	params theta;
	//Body Weight (kg)
	//cout<<"# Tablets per dose: "<<TABLETS<<endl;
	//Body Weight (kg)
	double BW;//11.1 + 2.8 * distribution(generator);//SHOULD BE RANDOM, AND DOSE SHOULD BE ADJUSTED ACCORDINGLY
	/*while(BW<5||BW>25){
		BW = 11.1 + 2.8 * distribution(generator);
	}*/
	
	if(TABLETS == 1){
		BW = 5 + 10 * distribution2(generator);
	}
	if(TABLETS == 2){
		BW = 15 + 10 * distribution2(generator);
	}
	if(TABLETS == 3){
		BW = 25 + 10 * distribution2(generator);
	}
	if(TABLETS == 4){
		BW = 35 + 10 * distribution2(generator);
	}
	//cout<<"Body weight: "<< BW <<endl;

	//Artemether constants
	theta.CLAM = 24.7 * pow(BW , 0.75) * exp(r1);//Write as function?
	theta.VCAM = 129;
	theta.VMAM = theta.VCAM;
	theta.kaAM = 0.27 * exp(r2);
	theta.k23AM = 5.86 * exp(r3);
	theta.CLmetAM = 419;

	//Lumefantrine constants
	theta.CLL = 0.84 * pow(BW , 0.52) * exp(r4);//Write as a function?
	theta.VCL = 59.9 * pow(BW , 0.35) * exp(r5);
	theta.VML = theta.VCL;
	theta.CLmetL = 4.8;
	theta.kaL = 0.54;
	theta.k23L = 3.7 * pow(10,-4) * exp(r6);

	//////////////////////////////////////////////////
	/* Generate PD parameters in the struct. */
	//////////////////////////////////////////////////

	paramsPD2 thetaPD2;
	//fix DHA values to AM values
	thetaPD2.kmax_AM = 0.189;
	thetaPD2.c50_AM = 3.3;

	thetaPD2.kmax_L = 0.165;
	thetaPD2.c50_L = 125.0;

	//NOTE: DLF effect turned off at present
	thetaPD2.kmax_DLF = 0.0 / 24;
	thetaPD2.c50_DLF= 280.0;

	thetaPD2.Gkmax_AM.push_back(0.18);
	thetaPD2.Gkmax_AM.push_back(0.18);
	thetaPD2.Gkmax_AM.push_back(0.11);
	thetaPD2.Gkmax_AM.push_back(0.11);
	thetaPD2.Gkmax_AM.push_back(0.08);

	thetaPD2.Gkmax_L.push_back(0.13);
	thetaPD2.Gkmax_L.push_back(0.13);
	thetaPD2.Gkmax_L.push_back(0.04);
	thetaPD2.Gkmax_L.push_back(0.04);
	thetaPD2.Gkmax_L.push_back(0.0);

	//cout<<"Length: "<<thetaPD2.Gkmax_L.size()<<" Check G1_LMF "<< thetaPD2.Gkmax_L[0]<<endl;

	//////////////////////////////////////////////////
	/* Declare vector to store parasitaemia, and give initial condition at t=0 */
	//////////////////////////////////////////////////
	vector<double> P(400,0.0);
	P[0]=0.1;//initial condition

	vector<double> Q(40,0.0);//vector to help assess Day0 & Day1 parasitaemia
	Q[0]=0.1;//initial condition

	vector<double> PC2q(40,0.0);
	double PCq;
	double PCsumq;
	double Pvq=0;

	double Scq = 1;
	double Smq = 1;
	double Svq = 1;
	int upperq, lowerq;
	double Qmax=0;//max, from untreated episode with same random numbers
	for(int k=0;k<20;k++){ //Loop for Q (untreated), not P

		Scq = 1/(1+pow((Q[k]/Pcs),kaC));
		if(k>3){
			lowerq = round((k+1-4)*pow(lambda , k-4+1))-1;//Minus 1, since in Mathematica Initial Condition is at k=1
			//cout<< (k-4)*pow(lambda , k-4) <<endl;
			upperq = k - 4 +1;
			//cout<<"Lower: "<<lower << " Upper: "<<upper <<endl;

			Pvq=0;
			for(int k2=lowerq;k2<upperq;k2++){
				Pvq += Q[k2];
			}
			Svq=1/(1+pow((Pvq/Pvs),kaV));
			}
		// General-adaptive immune response
		if(Q[k]>C){
			PCq=C;
		}
		else{
			PCq=Q[k];
		}

		PC2q[k]=PCq;
			
		if(k>3){
		
			PCsumq=0;
			for(int k1=0;k1<k-4+1;k1++){
				PCsumq+=PC2q[k1];
			}
			Smq = beta + (1-beta)/(1+pow((PCsumq/Pms),kaM));
		}

		Q[k+1] = R[k] * Scq * Smq * Svq * Q[k];
		if(Q[k+1] > Qmax){
			Qmax = Q[k+1];
		}
		if(Q[k+1]<pow(10,-5)){
			break;
		}
	}// end of Parasitaemia loop

	//////////////////////////////////////////////////
	/* Fever Threshold & Waiting time for treatment */
	//////////////////////////////////////////////////

	//predict the peak parasitaemia, based on (random) innate immunity factor
	double Ur = distribution3(generator);
	double thresh = pow(10,Ur + log10(Qmax) );
	//cout<<"Fever Threshold:"<< thresh <<endl;

	//Delay from fever breaking to 1st dose
	double meanLNw = 0.7;
	double sigmaLNw = 0.6575;

	double rand1w = meanLNw + sigmaLNw * distribution(generator);
	double rand2w = exp(rand1w);
	while(rand2w>10.0){//truncate dist
		rand1w = meanLNw + sigmaLNw * distribution(generator);
		rand2w = exp(rand1w);
	}
	int rand3w = round(rand2w * 24 *(1/double(dt) ));

	//cout<<"Delay to treatment (in Days): "<< rand2w <<endl;

	int uu1;
	int uu2;

	int AS = 0;
	int CF = 0;
	int fevk = 0;
	//If treatment fails, retreat with probability Pretreat
	double Pretreat = distribution2(generator);

	//////////////////////////////////////////////////
	/* Begin Runge-Kutta routine */
	//////////////////////////////////////////////////

	double k1L[3], k2L[3], k3L[3], k4L[3];
	double xL, yL, zL;
	//double k1AM[3], k2AM[3], k3AM[3], k4AM[3];
	int counterG;
	//double xAM, yAM, zAM;
	double killF = 0.0;
	double killFgam = 0.0;

	//Containers for the analytical solution for the Artemether PK ODE
	vector<double> store(3);
	vector<double> tempV(3);
	int cc;
	int summ = 0;

	for(int k = 0; k < TL - 1; k++){
		counterG = 0;
		DRUGL = 0.0;

		for(int k1=0;k1<24;k1++){
		if(k-DELTA==d24[k1] && d24[k1]>=0) //Remember entries include '-1's, for dose not taken or required
			counterG++; //More than one pill can be taken in one time step!
		}

		if(counterG>0){
			//DRUGAM = 20 * counterG;
			DRUGL = 120 * counterG;
			//cout<<"At increment "<< k <<" the dose is: "<< DRUGL <<" for LMF" <<endl;
		}
		/*else{
			DRUGL=0.0;
			}*/

		store[0] = 0.0;
		store[1] = 0.0;
		store[2] = 0.0;
		cc=0;

		if(k >= DELTA){

			for(int h=0;h<l1;h++){
				if(k-DELTA>=d11[h])
					cc++;
			}

			for(int j = 0; j < cc;j++){
				tempV  = ART(k*dt - DELTA*dt - d11[j]*dt,20*dm11[j],&theta);
				store[0] += tempV[0];
				store[1] += tempV[1];
				store[2] += tempV[2];
			}	
			/*if(i>=DELTA+392&&i<533)
				cout<<"Time: " <<i*dt<<" DELTA "<<DELTA*dt<<" cc equals "<<cc<<" "<<store[0]<<endl;*/
		}

		GUTAM[k] = store[0];
		CENTRALAM[k] = store[1];
		METABOLITEAM[k] = store[2];	

		//Convert to the desired units
		CENTRALAMx[k] = 1000 * (CENTRALAM[k]/theta.VCAM) ;
		METABOLITEAMx[k] = 1000 * (METABOLITEAM[k]/theta.VMAM) ;


		GUTL[k]+= DRUGL;	
		xL = GUTL[k];
		yL = CENTRALL[k];
		zL = METABOLITEL[k];

		k1L[0] = dt * f1L(xL, yL, zL, &theta);
		k1L[1] = dt * f2L(xL, yL, zL, &theta);
		k1L[2] = dt * f3L(xL, yL, zL, &theta);

		k2L[0] = dt * f1L(xL + 0.5 * k1L[0], yL + 0.5 * k1L[1], zL + 0.5 * k1L[2], &theta);
		k2L[1] = dt * f2L(xL + 0.5 * k1L[0], yL + 0.5 * k1L[1], zL + 0.5 * k1L[2], &theta);
		k2L[2] = dt * f3L(xL + 0.5 * k1L[0], yL + 0.5 * k1L[1], zL + 0.5 * k1L[2], &theta);

		k3L[0] = dt * f1L(xL + 0.5 * k2L[0], yL + 0.5 * k2L[1], zL + 0.5 * k2L[2], &theta);
		k3L[1] = dt * f2L(xL + 0.5 * k2L[0], yL + 0.5 * k2L[1], zL + 0.5 * k2L[2], &theta);
		k3L[2] = dt * f3L(xL + 0.5 * k2L[0], yL + 0.5 * k2L[1], zL + 0.5 * k2L[2], &theta);

		k4L[0] = dt * f1L(xL + k3L[0], yL + k3L[1], zL + k3L[2], &theta);
		k4L[1] = dt * f2L(xL + k3L[0], yL + k3L[1], zL + k3L[2], &theta);
		k4L[2] = dt * f3L(xL + k3L[0], yL + k3L[1], zL + k3L[2], &theta);

		GUTL[k+1] = GUTL[k] +  (k1L[0] + 2 * k2L[0] + 2 * k3L[0] + k4L[0])/6;
		CENTRALL[k+1] = CENTRALL[k] + (k1L[1] + 2 * k2L[1] + 2 * k3L[1] + k4L[1])/6;
		METABOLITEL[k+1] = METABOLITEL[k] + (k1L[2] + 2 * k2L[2] + 2 * k3L[2] + k4L[2])/6;

		//Convert to the desired units!
		CENTRALLx[k+1] = 1000 * (CENTRALL[k+1]/theta.VCL);
		METABOLITELx[k+1] = 1000 * (METABOLITEL[k+1]/theta.VML);

		//update 'kill factor'
		killF -= dt * ( thetaPD2.kmax_L * (CENTRALLx[k]/(CENTRALLx[k]+thetaPD2.c50_L)) /*+ 
		thetaPD.kmax_DLF * (METABOLITELx[k]/(METABOLITELx[k]+thetaPD.c50_DLF)) */+
		thetaPD2.kmax_AM * (CENTRALAMx[k]/(CENTRALAMx[k]+thetaPD2.c50_AM)) + 
		thetaPD2.kmax_AM * (METABOLITEAMx[k]/(METABOLITEAMx[k]+thetaPD2.c50_AM)) );

		fG1[0] = dt * gfPD(0,G[0][k],0,nuv[0],0,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2);
		fG2[0] = dt * gfPD(0,G[0][k] + 0.5 * fG1[0],0,nuv[0],0,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2);
		fG3[0] = dt * gfPD(0,G[0][k] + 0.5 * fG2[0],0,nuv[0],0,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2);
		fG4[0] = dt * gfPD(0,G[0][k] + fG3[0],0,nuv[0],0,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2);
		//RK4 step. Combine these four fns together
		G[0][k+1] = G[0][k] + (fG1[0] + 2 * fG2[0] + 2 * fG3[0] +fG4[0])/6;
		for(int a=1; a < L * 5; a++){
			uu1 = (a-1)/L;
			uu2 = a/L;
			fG1[a] = dt * (gfPD(G[a-1][k],G[a][k],nuv[uu1],nuv[uu2],uu2,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2)); //NOT a, CONVERT to stage!!
			fG2[a] = dt * (gfPD(G[a-1][k] + 0.5 * fG1[a-1],G[a][k] + 0.5 * fG1[a],nuv[uu1],nuv[uu2],uu2,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2));
			fG3[a] = dt * (gfPD(G[a-1][k] + 0.5 * fG2[a-1],G[a][k] + 0.5 * fG2[a],nuv[uu1],nuv[uu2],uu2,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2));
			fG4[a] = dt * (gfPD(G[a-1][k] + fG3[a-1],G[a][k] + fG3[a],nuv[uu1],nuv[uu2],uu2,CENTRALAMx[k],METABOLITEAMx[k],CENTRALLx[k],&thetaPD2));
			//RK4 step. Combine these four fns together
			G[a][k+1] = G[a][k] + (fG1[a] + 2 * fG2[a] + 2 * fG3[a] +fG4[a])/6;
		}
		//cout<<CENTRALLx[k] << " "<< CENTRALAMx[k]<< " "<< METABOLITEAMx[k]<<endl;

		//for(int k=0;k<399;k++){ 
		if(k%(48 * delta_t)==twoday-1 && AS==0){ //CHANGE TO 1/dt
			int k1 = ((k+1) * dt /48 )-1;
			//cout<<"Check k1: "<<k1<<endl;

			//Innate immune response
			Sc = 1/(1+pow((P[k1]/Pcs),kaC));

			//Effective var-specific adaptive immune response
			//First determine the limits on the sum?
			if(k1>3){
				lower = round((k1+1-4)*pow(lambda , k1-4+1))-1;//Minus 1, since in Mathematica Initial Condition is at k=1
				//cout<< (k-4)*pow(lambda , k-4) <<endl;
				upper = k1 - 4 +1;
				//cout<<"Lower: "<<lower << " Upper: "<<upper <<endl;

				Pv=0;
				for(int k2=lower;k2<upper;k2++){
					Pv += P[k2];
				}
				Sv=1/(1+pow((Pv/Pvs),kaV));
				//cout<<"k1 equals: "<<k1<<" Lower: "<<lower << " Upper: "<<upper <<" Pv: "<<Pv<<endl;
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

			P[k1+1] =  R[k1] * Sc * Sm * Sv * P[k1] * exp(killF);
			G[0][k+1] += alpha * R[k1] * Sc * Sm * Sv * P[k1] * exp(killF); //SHOLLD GROWTH RATE BE INCLUDED HERE?????
			//cout<<"On Day "<< 2*k1<<", " << P[k1]<<" "<<Sc<<" "<<Sm<<" "<<Sv<<" "<< G[0][k+1]<<endl;

			if(P[k1+1]>thresh && CF==0){
				DELTA = k + rand3w;
				fevk = k;
				//cout<<"Fever triggered! "<<k*dt/24<<", DELTA: "<<DELTA*dt/24<<endl;
				CF = 1;
			}

			//What about re-treatment?? Now f_last replaces D6 for final dose
			if(P[k1+1]>thresh && CF==1 && k> f_last + DELTA){
				if(Pretreat < prob_retreat){
					DELTA_old = DELTA; //store time of original treatment
					DELTA = k + rand3w;
					//cout<<"Woah, fever triggered again! "<<k*dt/24<<", DELTA: "<<DELTA*dt/24<<endl;
					CF = 2;
				}
				else{
					//cout<<"Fever triggered, but no treatment"<<endl;
					CF = 3;
				}
			}

			if(P[k1+1]<pow(10,-5)){
				AS=1;
			}

			killF = 0.0; //reset
		}//End of asexual parasitaemia loop

		if(AS==1){
			double bavuma = 0.0;
			for(int ii=0;ii<L;ii++){
				bavuma += G[4*L + ii][k];
			}

			if(AS>0 && bavuma < pow(10,-5)){
				//cout<< "Last mature gametocytes have gone" <<endl;
				break;
			}
		}


	}//end of 'dt' loop

	vector<double> GT1(TL,0.0);
	vector<double> GT2(TL,0.0);
	vector<double> GT3(TL,0.0);
	vector<double> GT4(TL,0.0);
	vector<double> GT5(TL,0.0);
	double P0, G0;//Asexual density & mature gametocytes at time of first dose

	for(int k=0;k<TL;k++){
		for(int i=0;i<L;i++){
			GT1[k] += G[i][k];
			GT2[k] += G[L + i][k];
			GT3[k] += G[2*L + i][k];
			GT4[k] += G[3*L + i][k];
			GT5[k] += G[4*L + i][k];
		}
	}
	G0 = GT5[DELTA];

	G.clear();

	///////////////////////////////////////////////////
	/*				AUC Calculation					*/
	//////////////////////////////////////////////////

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
	AUICvec[GG] += (1.0/NN)*auic;
} //q loop
cout<<"Mean AUIC for patient #"<<GG<<": "<<AUICvec[GG]<<endl;
}//End of patient loop

//Overall mean
double summ=0.0;
for(int j=0;j<659;j++){
	summ += (1.0/482)*AUICvec[j];
} 
cout<<"Overall, mean AUIC was"<<summ<<endl;
	///////////////////////////////////////////////////
	/*					Output 						*/
	//////////////////////////////////////////////////

	 	/*{
		  ofstream out1(output_File1);
		  if(!out1){
		    cerr <<"Failed to open output file"<< output_File1 << endl;
		    exit(1);
		  }*/
		  /*for(int f3=0;f3<30-1;f3++){
		    out1 << 2*f3 <<"\t"<<P[f3]<<"\t"<<Q[f3]<<"\t"<<
		    GT1[f3*twoday]<<"\t"<<
		    GT2[f3*twoday]<<"\t"<<
		    GT3[f3*twoday]<<"\t"<<
		    GT4[f3*twoday]<<"\t"<<
		    GT5[f3*twoday]<<endl;
		  }*/
		    //out1 << P0 <<"\t"<< G0 <<endl;
		//}

		/*{
		  ofstream out2(output_File2);
		  if(!out2){
		    cerr <<"Failed to open output file"<< output_File2 << endl;
		    exit(1);
		  }*/
		 
		// for(int f3 = 0;f3 < TL-1;f3+=5){
		//    out2<< dt*f3 <<"\t" << GT5[f3]<<"\t"<< P_infect[f3] <<endl;
		// }
		// //}

		// out1 << DELTA <<endl;
		// out1 << DELTA_old <<endl;
		// out1 << rand2w <<endl;
		// out1 << thresh <<endl;
		// for(int f3 = 0;f3 < 200;f3++){
		//     out1 << P[f3] <<endl;
		// }

		// for(int f3 = 0;f3 < TL-1;f3+=5){
		//     out3 << dt*f3 <<"\t"<<CENTRALAMx[f3]<<"\t"<<METABOLITEAMx[f3]<<"\t"<<
		//     CENTRALLx[f3]<<"\t"<<METABOLITELx[f3]<<endl;
		// }


  } //end of main
