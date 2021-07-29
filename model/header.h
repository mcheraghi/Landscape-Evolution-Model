#include<iostream>
#include<iomanip>
#include<fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <deque>
#include <ctime>
#include <cstdlib>
#include <list>
#include <sstream>
#include <unistd.h>
#include <chrono>
#include <random>
#include <thread>
#include <tuple>
using namespace std;
using namespace std::this_thread; // sleep_for, sleep_until
using namespace std::chrono; // nanoseconds, system_clock, seconds

#define M		121//exp8: 241(4mm), 121(8mm)			 exp 7:123 (8 mm)//244 (4 mm)
#define N          	239 //exp8: 477(4mm), 239(8mm)			 exp 7:233 (8 mm)//464 (4 mm)
#define dx 		8
#define dy 		dx
#define delta 		dx
#define mainSlope 	0.05
#define aveRain 	84.0 //(mm/h)


#define CFL 		0.9
#define Tmax            16.0
#define PI 		3.14159265359



struct objectives
{
	double a;
	double b;
	double c;
	double d;
};

class interface
{
	
	public:
	int runOption;
	void getRunOption(string addressRUN, string addressOPT);
};

class matrixOperations
{
	
	public:
	int transpose(double *A, double *B, int row, int col);
	int multiple(double *A, double *B, double *C, int Ar, int Ac, int Br, int Bc);
	int getMinor(double **src, double **dest, int row, int col, int order);
	double calcDeterminant( double **mat, int order);
	void matrixInversion(double *A, int order, double *Y);
	friend class functionSet;
	
};



class functionSet:public matrixOperations
{
	private:
		int *DIR, *INPUT, *drainTo, *drainFrom;
		double *CURVE, *R, *D, *alpha, *down, *up, *right, *left,*downBC, *upBC, *rightBC, *leftBC, *scanTime, *param,*Omega,*K;
		double *a, *b, *c, *d; 
		ifstream input1,input2;
		ofstream output1,output2, output3;
		
		int PHICalc(double *Q, double *SLOPE, double *PHI);
		int Thomas(double *,double *,double *,double *,double *,int);
		void d2s(double num, char s[], unsigned int width);
		
	public:
	
		
	double E, Dave, Kave, m, n, OmegaAve,noiseAve;//E:tectonic, Dave:Diffusion, K:convection, m:discharge expoenent,n:slope exponent,Tave:critical stream power (Omega =Tave*random number) , noiseAve:coefficient of deposition term	
	double t;
	int BCstep;
	int varD;
		double slopeLocal[8], dt, error;
		double epsilon, gamma; //epsilon: the value for filling process, gamma: the value for BC
		int allocate(double *X,double *Y,double *Z,double *Q,double *dZstar,double *Zstar,double *PHI);
		void readLand(double *X,double *Y,double *Z, string str1);
		
		int filling(double *, double *, double *);
		void smoothZ(double *);
		void direction(double *, double *, double * );
		void L_down(double *);
		void drainageQ(double *Q);
		int mainRiver(double *X, double *Y, double *Z, double *L, double *BV, double BVmin,  double *SLOPE, double *Q, double *CURVE,double *LL, double *PHI, string folder, string fileName);//the main path is determined base on the based value (BV) and its minimum value (BVmin)
		int calcDrainFrom(double *BV);
		void UpLenght(double *L);
		int statistics(double *X, double *Y, double *Z, double *Q, double *AREA, double *SLOPE, double *CURVE, double *PHI, double *L, double *LL, double percentage, string folder);
		void prepareFile(string folderName, string fileName);
		void drainageAREAandLL(double *AREA, double *LL, double *order); //The drainage area without considering the rainfall
		void writeExceedance(double *paramMain, double *paramTemp, int NumPoints, string folderName, string fileName); //it writes the exceedance probabilities of a specific parameter
		void Kmean(double *X, double *Y, double *center, int numpoints, int myK, string folderName, string fileName);
		void Ktest(double *X, double *Y, double *center, int maxK, string folderName, string fileName);
		void findMinMax(double*param, double *min, double *max);
		void normalize(double *X);
		double average(double *X);
		double RMS(double *X, int size);
		double correlation(double *X, double *Y, int size);
		double meanSquareError(double *X, double *Y, int size);
		int Qcomponents(double *X, double *Y, double *Q, double *Qx, double *Qy);

		void defineRain(string str1, int uni);
		void defineD(int varD);
		void defineRandomD(int varD);// 0 or 1 is put directly in the function
		void defineK(int varK);
		void defineOmega(int varOmega);

		void slopeCalc(double *Z, double *SLOPE);
		void curveCalc(double *SLOPE, double *CURVE);

		objectives multiObjective_Zerror_discharge(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string initialTime, string finalTime,string* filename);
		double objective_discharge(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string initialTime, string finalTime,string* filename);
		int temporalBCset(double *X, double *Y, double *Z,string* filename);	//sets the boundary according to the experimental data
		void filter33(double *param, double *filtered, int filtTime);
		double filternn(double *param, double *filtered, int n);//it filters param by an n by n cells
		double filterCELL(double *param, int i,int j, int n);//it filters param[i,j] by an n by n cells
		void filterNetwork(double *Z1, double *Z2, double *Y, double *Qopt);
		int dtCalc(double *Q, double *SLOPE);
		int RungeKutta(double *Z, double *Zstar, double *Q, double *SLOPE, double *dZstar, double *PHI, double *Y);
		int crankNicholson(double *Zstar, double *Y);	
		int boundarySet(double *Z, double *Y);
		int replaceError(double *Z, double *Zstar);
		int addNoise(double *Z, double dt); //in numeric file
		int writeDATA(double *X, double *Y, double *Z, double *X1, double *X2, double *X3, double *X4, double *X5, double *X6, double *X7, string folder, int step);
		int writeXYZ(double *X, double *Y, double *Z, string folder, int step);
       		int analytic(double *X, double *Y, double *Z);
		
		int RUN(double *X, double *Y, double *Z, double *Zstar, double *Q,  double *Qx,  double *Qy, double *SLOPE, double *CURVE, double *PHI, double *dZstar, double *AREA, double *LL,double *downL, string addressRUN, double initialTime, double finalTime, string initialTimeZ, string finalTimeZ);
		double objective1(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string initialTime, string finalTime,string* filename);
		double stepRun(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, double initialT, double finalT, string finalZ) ;
		double objective2(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar);
		int calibratePSO(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string addressOPT, string initialTime, string finalTime,string* filename);
		int calibrateMontCarlo(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI);
		int confidenceInterval(double q[], double percentage[], double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar,double T1, double T2, string initialTime, string finalTime);
		int confidenceInterval1(double q0[], double percentage[], double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar,double T1, double T2, string initialTime, string finalTime);
		void SetCursorPos(int XPos, int YPos);

};
		

