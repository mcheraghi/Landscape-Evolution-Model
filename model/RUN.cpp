#include "header.h"


int functionSet::RUN(double *X, double *Y, double *Z, double *Zstar, double *Q, double *Qx, double *Qy, double *SLOPE,double *CURVE, double *PHI, double *dZstar, double *AREA, double *LL, double *downL, string addressRUN, double initialTime, double finalTime, string initialTimeZ, string finalTimeZ)
{
	double *order;
	double random,rand1;
	double noise;
	default_random_engine generator;
	normal_distribution<double> distribution(0,0.503702461771);//The noise distribution
	std::cauchy_distribution<double> cauchy(0,0.11757900431008157);
  
	int i,j;
	order      = (double*)calloc(M*N , sizeof(double));
	//check the writedata function, it has been changed probably
	string makefolder = "mkdir " + addressRUN;
	system(makefolder.c_str());
	
	//now it is constant rainfall, when you change the setting delete this sentence
	double epsilon0 = 0.1;
	epsilon = epsilon0;
	t = initialTime;
	double step = (finalTime -  initialTime)/10.0; 
	double captureTime = t *2.0;
	int nn;
	for (nn = 0; nn < 7; nn++)
		if( initialTime == scanTime[nn])
			BCstep = nn;
	
	double BCtime =  scanTime[BCstep + 1] ;
		
	defineD(varD);
	defineK(0);
	defineOmega(1);
	readLand(X, Y, Z,initialTimeZ);///
	
	//filter33(Z, Z, 1);
	//writeDATA(X, Y, Z, R, D, addressRUN, 1000);
	replaceError(Zstar, Z);
	
	boundarySet(Zstar, Y);
	boundarySet(Z, Y);
	direction(Zstar, PHI, Y);
	drainageQ(Q);
	slopeCalc(Z,SLOPE);
	//we use the original data to calculate slope (without filling)
	double ERROR = 5;
	int num = 0;

	//drainageAREAandLL(AREA, LL);
	curveCalc(Z,CURVE);
	drainageAREAandLL(AREA, LL,order);
	writeDATA(X, Y,Z,Z,SLOPE,Q,AREA,LL, CURVE,downL, addressRUN, num);
	writeXYZ(X, Y, Z, addressRUN, num);

	
	while(t < finalTime )
	{
		boundarySet(Z, Y);
		
		dtCalc(Q, SLOPE);
		epsilon = epsilon0-(epsilon0 - epsilon0/10.0)/(finalTime - initialTime)*t;//pow(10,rndCoef)/2.0;//epsilon0-(epsilon0 - epsilon0/10.0)/(finalTime - initialTime)*t;	//is used in filling processes
		if(t + dt > captureTime)
			dt = captureTime - t;
		if(t + dt >  BCtime)
		{
			BCstep++; BCtime = scanTime[BCstep + 1];		
		}
		RungeKutta(Z, Zstar, Q, SLOPE, dZstar, PHI, Y);
		crankNicholson(Zstar, Y);
	

		replaceError(Z, Zstar);
		
		//addNoise(Z,dt);
		defineOmega(1);
		

		boundarySet(Z, Y);
		boundarySet(Zstar, Y);
		
		t+=dt;
		//ERROR = error;
		//cout <<"CHANGE (MM) =" << ERROR << "    " << "t =" << t << endl;
		cout << "t =" << t << endl;
		if(t == captureTime)
		{
			direction(Zstar, PHI, Y);//these two function are also in numeric PHIcalc
			drainageQ(Q);
			slopeCalc(Z,SLOPE);
			curveCalc(Z,CURVE); Qcomponents(X, Y, Q, Qx, Qy); 
			drainageAREAandLL(AREA, LL,order);		
			writeDATA(X, Y,Z,Z,SLOPE,Q,AREA,LL, CURVE, downL,addressRUN, ++num);
			writeXYZ(X, Y, Z, addressRUN,  num);
			captureTime*=2.0;
			
		}
		
		direction(Zstar, PHI, Y);//these two function are also in numeric PHIcalc
		drainageQ(Q);
		slopeCalc(Z,SLOPE);
			
	}
	readLand(X, Y, Zstar, finalTimeZ);
	

	ERROR = 0;
	cout <<"I AM HERE" << endl;
	ERROR = correlation(Zstar, Z, M*N);
	cout << "-----------ERROR = " << ERROR <<endl<<endl ;
	std::string command = "python PlotScript.py";
    	system(command.c_str());
	return 1;
}


