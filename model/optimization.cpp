#include "header.h"



double functionSet::stepRun(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, double initialT, double finalT, string finalZ) //Z is the initial Z
{

	epsilon = 1.0;
	double random,rand1;
	double noise;

	int i,j;

	//defineD(varD);
	defineRandomD(1);
	defineK(0);
	defineOmega(0);
	replaceError(Zstar, Z);
	direction(Zstar, PHI, Y);
	drainageQ(Q);
	slopeCalc(Z,SLOPE); //we use the original data to calculate slope (without filling)
	
	t = initialT;

	int nn;
	for (nn = 0; nn < 7; nn++)
		if( initialT == scanTime[nn])
			BCstep = nn;
	
	double BCtime =  scanTime[BCstep + 1] ;

	while(t < finalT)
	{
		//epsilon = 0.05-(0.05 - 1e-2)/finalT*t;	
		//cout << "    " << "t =" << t << endl;
		dtCalc(Q, SLOPE);
		epsilon = 0.05-(0.05 - 1e-2)/finalT*t; //rndCoef*pow(10,0.5)*dt
		if (dt < 0.5e-3)
			return 1/dt;
		if (dt + t > finalT)
			dt = finalT - t;

		if(t + dt >  BCtime)
		{
			BCstep++; BCtime = scanTime[BCstep + 1];		
		}
		RungeKutta(Z, Zstar, Q, SLOPE, dZstar, PHI, Y);
		crankNicholson(Zstar, Y);
		replaceError(Z, Zstar);

		defineRandomD(1);

		boundarySet(Z, Y);
		boundarySet(Zstar, Y);
		t+=dt;
		//cout << "    " << "t =" << t << endl;
		if (Kave == 0.0)
			continue;
		direction(Zstar, PHI, Y);
		drainageQ(Q);
		slopeCalc(Z,SLOPE);	
	}

	readLand(X, Y, Zstar,finalZ);
	
	
	return meanSquareError(Z, Zstar, M*N);

	
} 

objectives functionSet::multiObjective_Zerror_discharge(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string initialTime, string finalTime,string* filename)
{
	int i,j,k;
	double *Qs,*Ls;
	double rmsZ, rmsQ, rmsL;	
	objectives ERROR;
	ERROR.a = 0;
	ERROR.b = 0;
	ERROR.c = 0;
	ERROR.d = 0;
	readLand(X, Y, Z,initialTime);
	//MCH.smoothZ(Z);
	//cout<<"ok"<<endl;
	double time[] = {0.25,0.5,1.0,2.0,4.0,8.0};
	Qs = new double[6*M*N]; 
	Ls = new double[6*M*N];
	for (k = 1; k<6; k++)
	{
		rmsZ =  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,time[k-1], time[k], filename[k]); //0:uniform D
		drainageAREAandLL(SLOPE, PHI,dZstar); //because we don't need SLOPE and PHI at this step, they are used to save AREA and LL. AREA:SLOPE, LL:PHI, order:dZstar
		for(j = 0; j<N; j++)
		for(i = 0; i<M; i++)
		{
			Qs[k*M*N+j*M + i] = Q[j*M + i];	
			Ls[k*M*N+j*M + i] = PHI[j*M + i];	
		}
	}
	ERROR.a = rmsZ;
	output1.close();
	output1.open("Q_model.csv",ios::trunc);
	
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
	{
		for (k = 1; k<6; k++)
			output1 << Qs[k*M*N+j*M + i] <<',';
		output1 << endl;	
	}

//calculating Q rms
	output1.close();
	//Running the python script to calculate rms of discharges
	sleep_for(seconds(4));
	std::string command = "python Q_RMS.py";
    	system(command.c_str());
	sleep_for(seconds(4));
	input1.open("Q_RMS.dat");//
	input1>>rmsQ;
	ERROR.b = rmsQ;
	input1.close();
	remove("Q_RMS.dat");
	remove("Q_model.csv");


//calculating L rms
	output1.close();
	output1.open("L_model.csv",ios::trunc);
	
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
	{
		for (k = 1; k<6; k++)
			output1 << Ls[k*M*N+j*M + i] <<',';
		output1 << endl;	
	}
	output1.close();
	//Running the python script to calculate rms of discharges
	sleep_for(seconds(4));
	command = "python L_RMS.py";
    	system(command.c_str());
	sleep_for(seconds(4));
	input1.open("L_RMS.dat");//
	input1>>rmsL;
	ERROR.c = rmsL;
	input1.close();
	remove("L_RMS.dat");
	remove("L_model.csv");

	ERROR.d = 0;//rmsL*rmsQ*rmsZ;

	return ERROR;	
}




//Objective function is based on the exceedance probability of discharge at 8h 
double functionSet::objective_discharge(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string initialTime, string finalTime,string* filename)
{
	int i,j;
	readLand(X, Y, Z,initialTime);
	//MCH.smoothZ(Z);
	double ERROR = 0;
	double error1;
	//error1 =  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,0.25, 8.0, filename[5]); //0:uniform D
	
	//Writing the model Q data in the file
	output1.close();
	output1.open("model_discharge.csv",ios::trunc);
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
	{
		output1 << Q[j*M + i] <<endl;	
	}
	output1.close();
	//Running the python script to calculate rms of discharges
	sleep_for(seconds(4));
	std::string command = "python Q_RMS.py";
    	system(command.c_str());
	sleep_for(seconds(4));
	input1.open("Discharge_RMS.dat");//
	input1>>ERROR;
	input1.close();
	remove("Discharge_RMS.dat");
	remove("model_discharge.csv");
	return ERROR;
	
}


double functionSet::objective1(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string initialTime, string finalTime,string* filename)
{

	
	readLand(X, Y, Z,initialTime);
	//MCH.smoothZ(Z);
	double ERROR = 0;
	
	double error1 = 0;
	error1 =  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,0.25, 0.5, filename[1]); //0:uniform D
	output1 << "    "<< error1;
	//ERROR += error1;

	error1=  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,0.5, 1.0, filename[2]);
	output1 << "    "<< error1;
	//ERROR += error1;

	error1=  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,1.0, 2.0, filename[3]);
	output1 << "    "<< error1;
	//ERROR += error1;

	error1=  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,2.0, 4.0, filename[4]);
	output1 << "    "<< error1;
	//ERROR += error1;

	error1=  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,4.0, 8.0, filename[5]);
	output1 << "    "<< error1;
	ERROR += error1;

	//error1 =  stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,8.0, 16.0, filename[6]);
	//output1 << "    "<< error1;
	//ERROR += error1;
	
	//ERROR += stepRun(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar,0.25, 16.0, "/home/mohsen/LandScapeCode/input/exp8/16h00min_8mm.dat");//number "1" is for var D and "0" for fixed D
	
	return ERROR;
	
}


double functionSet::objective2(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar)
{
	int i,j;
	/*epsilon = 0.05;
	
	readLand(X, Y, Z,"/home/mohsen/LandScapeCode/input/exp8/16h00min_8mm.dat");
	//MCH.smoothZ(Z);
	boundarySet(Z, Y);
	 //It will be changed during the filling process, it is just the initial value
	replaceError(Zstar, Z);
	direction(Zstar, PHI, Y);
	drainageQ(Q);
	slopeCalc(Zstar,SLOPE);
	error = 0;
	double error1 = 0;
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
	{
		error1 = 0;
		if(pow(Q[j*M + i], m)*pow(SLOPE[j*M + i], n) > tetaC)
			error1 += E- delta/dx*(pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) - tetaC);			
		error += error1*error1;
	}
	
	return sqrt(error/(M*N));*/

	
	epsilon = 0.05;
	defineD(varD);
	readLand(X, Y, Z,"/home/mohsen/LandScapeCode/input/exp8/16h00min_8mm.dat");
	//MCH.smoothZ(Z);
	boundarySet(Z, Y);
	 //It will be changed during the filling process, it is just the initial value
	replaceError(Zstar, Z);
	direction(Zstar, PHI, Y);
	drainageQ(Q);
	slopeCalc(Zstar,SLOPE);
	
	error = 0;
	double error1 = 0;
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
	{
		error1 = 0;
		error1 += D[j*M + i]*(Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] - 4*Z[j*M + i])/(dx*dx) + (D[j*M + i+1] - D[j*M + i])* (Z[j*M + i + 1] - Z[j*M + i])/(dx*dx)
		+ (D[(j + 1)*M + i] - D[j*M + i])* (Z[(j + 1)*M + i] - Z[j*M + i])/(dy*dy);

		if(pow(Q[j*M + i], m)*pow(SLOPE[j*M + i], n) > OmegaAve)
			error1 += E - delta/dx*K[j*M + i]*(pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) - Omega[j*M + i]);
					
		error += error1*error1;
	}
	
	return sqrt(error/(M*N));
	
	
}



int functionSet::calibratePSO(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar, string addressOPT, string initialTime, string finalTime,string* filename)
{
	
	int i,j;
	remove(addressOPT.c_str());
	output1.open(addressOPT.c_str(),ios::app);
	objectives objFun;
	int nvar = 7;


	double min[nvar];
	min[0] = 16000 ;//D, 
	min[1] = 0.17;//K, 
	min[2] = 10;//Tc
	min[3] = 0;//gamma, 
	min[4] = 0.5;//m, 
	min[5] = 1.0;//n, 
	min[6] = 0.0;//n, 
	
	double minOLD[nvar];
	double maxOLD[nvar];

	double max[nvar];
	max[0] = 21000;//D, 
	max[1] =  0.23;//K, 
	max[2] =10;//Tc 
	max[3] = 0;//gamma, 
	max[4] = 0.5;//m, 
	max[5] = 1.0;//n, 
	max[6] = 100.0;//randCoef, 
	E = 0; //E
	int pop = 50;
	
	

	double *Point, *obj, *selfOptPlace, *vel, *selfOpt;

	Point  = (double*)calloc(pop*nvar , sizeof(double));
	selfOptPlace  = (double*)calloc(pop*nvar , sizeof(double));
	vel  =     (double*)calloc(pop*nvar , sizeof(double));
	obj  =     (double*)calloc(pop , sizeof(double));
	selfOpt =  (double*)calloc(pop , sizeof(double));
	
	double optGeneral = 	1e6;

	double optGeneralPlace[] ={1000,1000,0, 0, 0.5, 1.0,100}; 
	
	
	int mm, nn;
	/*string fileNAme= "/home/mohsen/LandScapeCode/2D/Results/exp8/OPTIMIZATION/mn0.6/finalW.dat";
	input1.open(fileNAme.c_str());
	double numb;

	for(mm = 0; mm<pop; mm++) 
	{input1 >> numb;//initialization of data
		for(nn = 0; nn < 5; nn++)
		{
			input1 >> Point[mm*5 + nn];
			selfOptPlace[mm*5 + nn] = Point[mm*5 + nn];
			
		}
	input1 >> obj[mm]; 
			selfOpt[mm] =obj[mm];
	
	}
	for(mm = 0; mm<pop; mm++)
	{
		
	 << 0 <<"    "<<mm;	
		output1 <<mm << " ";	
		for(nn = 0; nn < 5; nn++)
			{
				output1 <<" "<< Point[mm*5 + nn];
				cout <<" "<< Point[mm*5 + nn];
			}	
		output1 << " "<<obj[mm]<<endl;
		cout << " "<<obj[mm]<<endl;	
	}*/
	double random;
	for(mm = 0; mm<pop; mm++) 
	{
		//random = (double)rand()/(double)RAND_MAX;
		for(nn = 0; nn < nvar; nn++)
		{
			random = (double)rand()/(double)RAND_MAX;
			Point[mm*nvar + nn] = min[nn] + random*(max[nn] - min[nn]) ; 	//D
			vel[mm*nvar + nn] = 0;
			selfOptPlace[mm*nvar + nn] = Point[mm*nvar + nn];
			selfOpt[mm] = 1e15;
		}
		vel[mm*nvar + 1] = 0;
		selfOptPlace[mm*nvar + 1] = Point[mm*nvar + 1];

		
	//Point[mm*nvar + 1] = 0.5*(double)rand()/(double)RAND_MAX*68.0601508637*pow(306179.135144142, -Point[mm*nvar + 2]);
	}
	

	for(nn = 0; nn < nvar; nn++)
		output1 << min[nn] <<"     "<<max[nn]<< endl;	
			
	
	int iter = 0;
	while ( iter < 100)
	{
		for(mm = 0; mm<pop; mm++)
		{
			
			cout << iter <<"    "<<mm;
			
			for(nn = 0; nn < nvar; nn++)			
				cout <<" "<< Point[mm*nvar + nn];
			
			output1 <<mm << " ";	
			for(nn = 0; nn < nvar; nn++)	
				output1 <<" "<< Point[mm*nvar + nn];

	
			Dave = Point[mm*nvar + 0];
			Kave = Point[mm*nvar + 1];
			OmegaAve = Point[mm*nvar + 2];
			gamma = Point[mm*nvar + 3];
			m = Point[mm*nvar + 4];
			n = Point[mm*nvar + 5];
			noiseAve = Point[mm*nvar + 6];
			E = 0;
			//obj[mm] = objective_discharge(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar, initialTime,finalTime,filename );
			//obj[mm] = objective2(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar);
			objFun = multiObjective_Zerror_discharge(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar, initialTime,finalTime,filename ); obj[mm] = objFun.a*objFun.b;
	
			if(obj[mm]<optGeneral)
			{
				optGeneral = obj[mm];
				for(nn = 0; nn < nvar; nn++)
					optGeneralPlace[nn] = Point[mm*nvar + nn];		
			} 

			if(obj[mm]<selfOpt[mm])
			{
				selfOpt[mm] = obj[mm];
				for(nn = 0; nn < nvar; nn++)
					selfOptPlace[mm*nvar + nn] = Point[mm*nvar + nn];		
			} 
			
			cout << " "<<objFun.a<<" "<<objFun.b<<" "<<obj[mm]<<endl;		
			output1 << " "<<objFun.a<<" "<<objFun.b<<" "<<obj[mm]<<endl;
				
			
			
		}

		for(nn = 0; nn < nvar; nn++)
		{
			for(mm = 0; mm<pop; mm++)
			{
				vel[mm*nvar + nn] = vel[mm*nvar + nn] + 2*(double)rand()/(double)RAND_MAX*1*(selfOptPlace[mm*nvar + nn] - Point[mm*nvar + nn]) + 2*(double)rand()/(double)RAND_MAX*1*(optGeneralPlace[nn] - Point[mm*nvar + nn]) ;
				Point[mm*nvar + nn] += vel[mm*nvar + nn];
		
				if (Point[mm*nvar + nn] > max[nn])
				{
					Point[mm*nvar + nn] = min[nn] + (double)rand()/(double)RAND_MAX*(max[nn] - min[nn]) ; //min[nn] + (double)rand()/(double)RAND_MAX*(max[nn] - min[nn]); //DO NOT use the limits point = max or min!!!!!
					//vel[mm*5 + nn]= 0;
				}
				if (Point[mm*nvar + nn] < min[nn])
				{
					Point[mm*nvar + nn] = min[nn] + (double)rand()/(double)RAND_MAX*(max[nn] - min[nn]);//min[nn] + (double)rand()/(double)RAND_MAX*(max[nn] - min[nn]); //DO NOT use the limits point = max or min!!!!!
					//vel[mm*5 + nn]= 0;	
				}
				
				
			}
		}
		
		
			
		
		cout <<endl<<endl<<endl<<iter << " "<< optGeneral << endl;
		output1 << iter<< " "<< optGeneral << endl;	
		iter++;
		for(nn = 0; nn < nvar; nn++)
		{
			output1 << optGeneralPlace[nn] << endl;	
			cout << optGeneralPlace[nn] << endl;
		}
		
		/*for(mm = 0; mm<pop; mm++)//remove this loop ;)
		{
			D = Point[mm*5 + 0];
			K = Point[mm*5 + 1];
			m = Point[mm*5 + 2];
			n = Point[mm*5 + 3];
			tetaC = Point[mm*5 + 4];
			E = 0;
			obj[mm] = objective1(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar, initialTime,finalTime );
			cout << iter <<"    "<<mm;	
			output1 <<mm << " ";	
			for(nn = 0; nn < 5; nn++)
				{
					output1 <<" "<< Point[mm*5 + nn];
					cout <<" "<< Point[mm*5 + nn];
				}	
			output1 << " "<<obj[mm]<<endl;
			cout << " "<<obj[mm]<<endl;
			
		}*/
	}
	
	/*int pop = 100.0;
	int iter = 0;
	remove(addressOPT.c_str());
	output1.open(addressOPT.c_str(),ios::app);
	double *Point, *obj, *selfOptPlace, *vel, *selfOpt;

	Point  = (double*)calloc(pop*5 , sizeof(double));
	selfOptPlace  = (double*)calloc(pop*5 , sizeof(double));
	vel  =     (double*)calloc(pop*5 , sizeof(double));
	obj  =     (double*)calloc(pop , sizeof(double));
	selfOpt =  (double*)calloc(pop , sizeof(double));
	
	double optGeneral1 = 1e15;
	double optGeneralPlace[5];
	double optGeneral2 = 1e15;
	optGeneralPlace[0] = 1074.75;
	optGeneralPlace[1] = 0.332554;
	optGeneralPlace[2] = 0.5;
	optGeneralPlace[3] = 1.0;
	optGeneralPlace[4] = 0.963411;
	

	while(iter < 100)
	{
		double max[5];
		max[0] = 2000.0;//D
		max[1] = 0.5;//K, 
		max[2] = 0.5;//m
		max[3] = 1.0;//n
		max[4] = optGeneralPlace[4];//Tc
	
		double min[5];
		min[0] = 500;//D
		min[1] = 0.2;//K,  
		min[2] = 0.5;//m
		min[3] = 1.0;//n
		min[4] = optGeneralPlace[4];//Tc

		
	
		E = 0;
		int mm, nn;
		for(mm = 0; mm<pop; mm++)
		for(nn = 0; nn < 5; nn++)
		{
			Point[mm*5 + nn] = min[nn] + (double)rand()/(double)RAND_MAX*(max[nn] - min[nn]) ; 	//D
			selfOptPlace[mm*5 + nn] = Point[mm*5 + nn];
			selfOpt[mm] = 1e15;	
		}
		
		for(mm = 0; mm<pop; mm++)
		{
			D = Point[mm*5 + 0];
			K = Point[mm*5 + 1];
			m = Point[mm*5 + 2];
			n = Point[mm*5 + 3];
			tetaC = Point[mm*5 + 4];
			E = 0;
			obj[mm] = objective1(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar );
	
			if(obj[mm]<optGeneral1)
			{
				optGeneral1 = obj[mm];
				for(nn = 0; nn < 5; nn++)
					optGeneralPlace[nn] = Point[mm*5 + nn];		
			} 

			if(obj[mm]<selfOpt[mm])
			{
				selfOpt[mm] = obj[mm];
				for(nn = 0; nn < 5; nn++)
					selfOptPlace[mm*5 + nn] = Point[mm*5 + nn];		
			} 

		
			cout << iter <<"    "<<mm;	
			output1 <<mm << " ";	
			for(nn = 0; nn < 5; nn++)
				{
					output1 <<" "<< Point[mm*5 + nn];
					cout <<" "<< Point[mm*5 + nn];
				}	
			output1 << " "<<obj[mm]<<endl;
			cout << " "<<obj[mm]<<endl;
			
		}


		cout <<endl<<endl<<endl<<iter << " "<< optGeneral1 << endl;
		output1 << iter<< " "<< optGeneral1 << endl;	

		for(nn = 0; nn < 5; nn++)
		{
			output1 << optGeneralPlace[nn] << endl;	
			cout << optGeneralPlace[nn] << endl;
		}
	



//at 20h
		optGeneral2 = 1e15;
		max[0] = optGeneralPlace[0];//D
		max[1] = optGeneralPlace[1];//K, 
		max[2] = 0.5;//m
		max[3] = 1.0;//n
		max[4] = 2*optGeneralPlace[4];//Tc
	
		min[0] = optGeneralPlace[0];//D
		min[1] = optGeneralPlace[1];//K,  
		min[2] = 0.5;//m
		min[3] = 1.0;//n
		min[4] = 0.5*optGeneralPlace[4];//Tc

		for(mm = 0; mm<pop; mm++)
		for(nn = 0; nn < 5; nn++)
		{
			Point[mm*5 + nn] = min[nn] + (double)rand()/(double)RAND_MAX*(max[nn] - min[nn]) ; 	//D
			selfOptPlace[mm*5 + nn] = Point[mm*5 + nn];
			selfOpt[mm] = 1e15;	
		}
		
		for(mm = 0; mm<pop; mm++)
		{
			D = Point[mm*5 + 0];
			K = Point[mm*5 + 1];
			m = Point[mm*5 + 2];
			n = Point[mm*5 + 3];
			tetaC = Point[mm*5 + 4];
			E = 0;
			obj[mm] = objective2(X, Y, Z, Zstar, Q, SLOPE,PHI, dZstar );
	
			if(obj[mm]<optGeneral2)
			{
				optGeneral2 = obj[mm];
				for(nn = 0; nn < 5; nn++)
					optGeneralPlace[nn] = Point[mm*5 + nn];		
			} 

			if(obj[mm]<selfOpt[mm])
			{
				selfOpt[mm] = obj[mm];
				for(nn = 0; nn < 5; nn++)
					selfOptPlace[mm*5 + nn] = Point[mm*5 + nn];		
			} 

		
			
			for(nn = 0; nn < 5; nn++)
				{
					cout <<" "<< Point[mm*5 + nn];
				}	
			
			cout <<endl;
		}
		cout <<endl<<endl<<endl<<iter << " "<< optGeneral2 << endl;	
		for(nn = 0; nn < 5; nn++)
		{
			cout << optGeneralPlace[nn] << endl;
		}
		iter++;
	}*/

	return 1;
}



void functionSet::SetCursorPos(int XPos, int YPos)
{
 printf("\033[%d;%dH", YPos+1, XPos+1);
}

int functionSet::confidenceInterval(double q[], double percentage[], double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar,double T1, double T2, string initialTime, string finalTime) //This function uses the value of z in the model to calculate the confidence interval using the method described in the uncertainty book
{
	int i,j;
	double *Z1, *Z2, *sens, *sensT;
	double a, sigma;
	Z1 	= new double[M*N];
	Z2 	= new double[M*N];
	sens 	= new double[M*N*3];
	sensT 	= new double[M*N*3];

	Dave 	= q[0];
	Kave 	= q[1];
	OmegaAve 	= q[2];
	m = 0.5;
	n = 1.0;
	
	readLand(X, Y, Z,initialTime);
	a = stepRun(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1, T2, finalTime);
	replaceError(Z1, Z);
	readLand(X, Y, Z,finalTime);
	sigma  = meanSquareError(Z1,Z,M*N/2)*sqrt(M*N/2)/sqrt(M*N/2 - 3);
	int nn, mm;
	cout << "D = 	" << q[0]<<"	"<<"	sigma = " << sigma << "		t = 1.9695 (95%)" <<endl;
	cout << "K = 	" << q[1]<<"	"<<"	sigma = " << sigma << "   	t = 1.9695 (95%)" <<endl;
	cout << "Tc = 	" << q[2]<<"	"<<"	sigma = " << sigma << "		t = 1.9695 (95%)" <<endl;
	cout  <<"dq1 %"<<"	dq2 %"<<"	dq3 %"<<"	sigmaD"<<"	sigmaK "<< "	sigmaTc    " << endl;

	double range = 0.05;
	double dq ;
	for(i = 0; i < 125; i++)
	{
		dq = -range+(double)rand()/(double)RAND_MAX*2*range;
		cout  << dq*100;
		Dave 	= q[0]*(1 + dq);
		Kave 	= q[1];
		OmegaAve 	= q[2];
		readLand(X, Y, Z,initialTime);
		a = stepRun(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1, T2, finalTime);
		replaceError(Z2, Z);

		for (mm = 0; mm < M*N/2; mm++)
		{
			sens[0*M*N + mm] = (Z2[mm] - Z1[mm])/(q[0]*dq);	
		}
		
		dq = -range+(double)rand()/(double)RAND_MAX*2*range;
		cout  << "    "<< dq*100;
		Dave 	= q[0];
		Kave 	= q[1]*(1 +dq);
		OmegaAve = q[2];
		readLand(X, Y, Z,initialTime);
		a = stepRun(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1, T2, finalTime);
		replaceError(Z2, Z);

		for (mm = 0; mm < M*N/2; mm++)
		{
			sens[1*M*N + mm] = (Z2[mm] - Z1[mm])/(q[1]*dq);	
		}

		dq = -range+(double)rand()/(double)RAND_MAX*2*range;
		cout  << "    "<< dq*100;
		Dave 	= q[0];
		Kave 	= q[1];
		OmegaAve 	= q[2]*(1 +dq);
		readLand(X, Y, Z,initialTime);
		a = stepRun(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1, T2, finalTime);
		replaceError(Z2, Z);

		for (mm = 0; mm < M*N/2; mm++)
		{
			sens[2*M*N + mm] = (Z2[mm] - Z1[mm])/(q[2]*dq);	
		}


		transpose(sens, sensT, M*N/2, 3);
	
		multiple(sensT, sens, Z, 3, M*N/2, M*N/2, 3);
		//for (int i = 0; i<3;i++)
		//	cout << "matrix:"<< Z[0*3 + i] <<"  "<<Z[1*3 + i] <<"  "<<Z[2*3 + i] <<endl;
		matrixInversion(Z, 3, SLOPE);
		cout <<"	"<<1.9695*sigma*sqrt(SLOPE[0*3 + 0]) <<"		"<<1.9695*sigma*sqrt(SLOPE[1*3 + 1]) <<"		"<< 1.9695*sigma*sqrt(SLOPE[2*3 + 2]) << endl;
	
	}
	output1.close();
	delete [] Z1;
	delete [] Z2;
	delete [] sens;
	delete [] sensT;
	
	
	return (0);

}


int functionSet::confidenceInterval1(double q0[], double percentage[], double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI, double *dZstar,double T1, double T2, string initialTime, string finalTime) //This function uses the objective function and second order term  to calculate the confidence interval (the method described by Prof.Davison)
{
	int Samp_Pop = 200;
	int i,j;
	double *Z1, *Z2, *sens, *sensT, *q,*sigma;
	double a, sigma0, sigmaBAR;
	Z1 	= new double[M*N];
	Z2 	= new double[M*N];
	q = new double[3*Samp_Pop];
	sigma = new double[Samp_Pop];
	sens 	= new double[Samp_Pop*6];
	sensT 	= new double[Samp_Pop*6];

	string str1 = "ConfidenceDATA_confid6percent.csv";
	output1.open(str1.c_str(),ios::app);

	Dave 	= q0[0];
	Kave 	= q0[1];
	OmegaAve 	= q0[2];
	m = 0.5;
	n = 1.0;
	
	readLand(X, Y, Z,initialTime);
	a = stepRun(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1, T2, finalTime);
	replaceError(Z1, Z);
	readLand(X, Y, Z,finalTime);
	sigma0  = meanSquareError(Z1,Z,M*N)*sqrt(M*N)/sqrt(M*N - 3);
	int nn, mm;
	cout << "D = 	" << q0[0]<<"	"<<"	sigma = " << sigma0 << "		t = 1.9695 (95%)" <<endl;
	cout << "K = 	" << q0[1]<<"	"<<"	sigma = " << sigma0 << "   	t = 1.9695 (95%)" <<endl;
	cout << "Tc = 	" << q0[2]<<"	"<<"	sigma = " << sigma0 << "		t = 1.9695 (95%)" <<endl;


	double range = 0.06;
	double dq[3] ;
	for(i = 0; i < Samp_Pop; i++)
	{
		dq[0] = -range+(double)rand()/(double)RAND_MAX*2*range;
		dq[1] = -range+(double)rand()/(double)RAND_MAX*2*range;
		dq[2] = -range+(double)rand()/(double)RAND_MAX*2*range;

		for (j = 0;j<3;j++)
			if (dq[j] ==0)
			{
				dq[j] = 1e-6;
			}
		q[0*Samp_Pop+i] = q0[0]*(1 + dq[0]);
		q[1*Samp_Pop+i] = q0[1]*(1 + dq[1]);
		q[2*Samp_Pop+i] = q0[2]*(1 + dq[2]);
		
		Dave 	= q[0*Samp_Pop+i] ;
		Kave 	= q[1*Samp_Pop+i] ;
		OmegaAve 	= q[2*Samp_Pop+i] ;

		readLand(X, Y, Z,initialTime);
		a = stepRun(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1, T2, finalTime);
		replaceError(Z1, Z);
		readLand(X, Y, Z,finalTime);
		sigma[i]  = meanSquareError(Z1,Z,M*N)*sqrt(M*N)/sqrt(M*N - 3);
		output1 <<q[0*Samp_Pop+i]<<"	"<<q[1*Samp_Pop+i]<<"	"<<q[2*Samp_Pop+i]<<"	"<<sigma[i]<<endl;
		cout <<q[0*Samp_Pop+i]<<"	"<<q[1*Samp_Pop+i]<<"	"<<q[2*Samp_Pop+i]<<"	"<<sigma[i]<<endl;
	}
	/*output1.close();  
	sigma0 = 10.4756;
	str1 = "/home/mohsen/LandScapeCode/2D/mainCode/ConfidenceDATA.csv";
	input1.open(str1.c_str());
	if (!input1)
	{
		cout << "The file can't be opened" << endl; 
		cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm" << endl;
		cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm" << endl;
		cout << "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMm" << endl;
		cout << str1<<endl;	
	}
	else
	{ 
		for(i = 0; i < Samp_Pop; i++)	
		{
		
			input1>>q[0*Samp_Pop+i];
			input1>>q[1*Samp_Pop+i];
			input1>>q[2*Samp_Pop+i];
			input1>>sigma[i];
		}
		input1.close();  
	}*/

	sigmaBAR = 0;
	for(i = 0; i < Samp_Pop; i++)
	{
		sigmaBAR += (sigma[i]-sigma0)*(sigma[i]-sigma0);
	}

	sigmaBAR = sqrt(sigmaBAR)/sqrt(Samp_Pop - 1);
	cout << "sigmaBAR="<<sigmaBAR<<endl;
	for(i = 0; i < Samp_Pop-1; i++)
	{
		
		sens[0*Samp_Pop+ i] =0.5*( (sigma[i]-sigma0)/(q[0*Samp_Pop+i]-q0[0]) - (sigma[i+1]-sigma0)/(q[0*Samp_Pop+i+1]-q0[0]))/(q[0*Samp_Pop+i]-q[0*Samp_Pop+i+1]);
		sens[1*Samp_Pop+ i] =0.5*( (sigma[i]-sigma0)/(q[1*Samp_Pop+i]-q0[1]) - (sigma[i+1]-sigma0)/(q[1*Samp_Pop+i+1]-q0[1]))/(q[1*Samp_Pop+i]-q[1*Samp_Pop+i+1]);
		sens[2*Samp_Pop+ i] =0.5*( (sigma[i]-sigma0)/(q[2*Samp_Pop+i]-q0[2]) - (sigma[i+1]-sigma0)/(q[2*Samp_Pop+i+1]-q0[2]))/(q[2*Samp_Pop+i]-q[2*Samp_Pop+i+1]);

		sens[3*Samp_Pop+ i] =( (sigma[i]-sigma0)/(q[0*Samp_Pop+i]-q0[0]) - (sigma[i+1]-sigma0)/(q[0*Samp_Pop+i+1]-q0[0]))/(q[1*Samp_Pop+i]-q[1*Samp_Pop+i+1]);
		sens[4*Samp_Pop+ i] =( (sigma[i]-sigma0)/(q[0*Samp_Pop+i]-q0[0]) - (sigma[i+1]-sigma0)/(q[0*Samp_Pop+i+1]-q0[0]))/(q[2*Samp_Pop+i]-q[2*Samp_Pop+i+1]);
		sens[5*Samp_Pop+ i] =( (sigma[i]-sigma0)/(q[1*Samp_Pop+i]-q0[1]) - (sigma[i+1]-sigma0)/(q[1*Samp_Pop+i+1]-q0[1]))/(q[2*Samp_Pop+i]-q[2*Samp_Pop+i+1]);
		cout <<i;
		for(j = 0;j<6;j++)
			cout<<"	"<<sens[j*Samp_Pop+ i];	
		cout <<endl;	
	}
	i = Samp_Pop - 1;
	sens[0*Samp_Pop+ i] =0.5*( (sigma[i]-sigma0)/(q[0*Samp_Pop+i]-q0[0]) - (sigma[i-1]-sigma0)/(q[0*Samp_Pop+i-1]-q0[0]))/(q[0*Samp_Pop+i]-q[0*Samp_Pop+i-1]);
	sens[1*Samp_Pop+ i] =0.5*( (sigma[i]-sigma0)/(q[1*Samp_Pop+i]-q0[1]) - (sigma[i-1]-sigma0)/(q[1*Samp_Pop+i-1]-q0[1]))/(q[1*Samp_Pop+i]-q[1*Samp_Pop+i-1]);
	sens[2*Samp_Pop+ i] =0.5*( (sigma[i]-sigma0)/(q[2*Samp_Pop+i]-q0[2]) - (sigma[i-1]-sigma0)/(q[2*Samp_Pop+i-1]-q0[2]))/(q[2*Samp_Pop+i]-q[2*Samp_Pop+i-1]);

	sens[3*Samp_Pop+ i] =( (sigma[i]-sigma0)/(q[0*Samp_Pop+i]-q0[0]) - (sigma[i-1]-sigma0)/(q[0*Samp_Pop+i-1]-q0[0]))/(q[1*Samp_Pop+i]-q[1*Samp_Pop+i-1]);
	sens[4*Samp_Pop+ i] =( (sigma[i]-sigma0)/(q[0*Samp_Pop+i]-q0[0]) - (sigma[i-1]-sigma0)/(q[0*Samp_Pop+i-1]-q0[0]))/(q[2*Samp_Pop+i]-q[2*Samp_Pop+i-1]);
	sens[5*Samp_Pop+ i] =( (sigma[i]-sigma0)/(q[1*Samp_Pop+i]-q0[1]) - (sigma[i-1]-sigma0)/(q[1*Samp_Pop+i-1]-q0[1]))/(q[2*Samp_Pop+i]-q[2*Samp_Pop+i-1]);
	
	transpose(sens, sensT, Samp_Pop, 6);
	
	multiple(sensT, sens, Z, 6, Samp_Pop, Samp_Pop, 6);
		//for (int i = 0; i<3;i++)
		//	cout << "matrix:"<< Z[0*3 + i] <<"  "<<Z[1*3 + i] <<"  "<<Z[2*3 + i] <<endl;
	matrixInversion(Z, 6, SLOPE);
	cout <<"	"<<1.9695*sigmaBAR*sqrt(SLOPE[0*6 + 0]) <<"		"<<1.9695*sigmaBAR*sqrt(SLOPE[1*6 + 1]) <<"		"<< 1.9695*sigmaBAR*sqrt(SLOPE[2*6 + 2]) << endl;
	output1 <<"	"<<1.9695*sigmaBAR*sqrt(SLOPE[0*6 + 0]) <<"		"<<1.9695*sigmaBAR*sqrt(SLOPE[1*6 + 1]) <<"		"<< 1.9695*sigmaBAR*sqrt(SLOPE[2*6 + 2]) << endl;
	
	
	delete [] Z1;
	delete [] Z2;
	delete [] sens;
	delete [] sensT;
	delete [] q;
	return (0);

}



/*
int functionSet::calibrateMontCarlo(double *X, double *Y, double *Z, double *Zstar, double *Q, double *SLOPE, double *PHI)
{
	

	readLand(X, Y, Z,"/home/mohsen/LandScapeCode/input/XYZ15%80mmh_4mmGrd_2h.txt");
	boundarySet(Z);
	epsilon = 0.05;   //the value critical for filling
	filling(Z, Zstar);
	direction(Z, Zstar);
	drainageQ(Q);
	slopeCalc(Z,SLOPE);
	
	readLand(X, Y, Zstar, "/home/mohsen/LandScapeCode/input/XYZ15%80mmh_4mmGrd_6h.txt");
	boundarySet(Zstar);
	filling(Zstar, PHI);
	direction(Zstar, PHI);
			
	
	
	remove("/home/mohsen/LandScapeCode/Results/calibrate/PSO/CalibratePSO2MILpoints.csv");
	output1.open("/home/mohsen/LandScapeCode/Results/calibrate/PSO/CalibratePSO2MILpoints.csv",ios::app);
	
	double max[5];
	max[0] = 50.0;
	max[1] = 50.0;
	max[2] = 1;
	max[3] = 1;
	max[4] = 50.0;
	
	
	double optGeneral = 1e15;
	double optGeneralPlace[5];
	int nn;
	int iter = 0;
	double obj;
	while ( iter < 2000000)
	{
		
			D = (double)rand()/(double)RAND_MAX*max[0]; 	//D;
			K = (double)rand()/(double)RAND_MAX*max[1];	//K;
			m = 0.5;					//m;
			n = 1.0;						//n
			tetaC = (double)rand()/(double)RAND_MAX*max[4];	//tetaC
			E = 0;
			obj = objective(Z,Zstar, Q, SLOPE);
	
			if(obj<optGeneral)
			{
				optGeneral = obj;
				optGeneralPlace[0] = D;	
				optGeneralPlace[1] = K;	
				optGeneralPlace[2] = m;		
				optGeneralPlace[3] = n;	
				optGeneralPlace[4] = tetaC;		
			} 

			

	
		if(iter % 10000 == 0)
		{
			cout << iter << "        "<< optGeneral << endl;
			output1 << iter<< "        "<< optGeneral << endl;	
			
			for(nn = 0; nn < 5; nn++)
				output1 << optGeneralPlace[nn] << endl;	
		}
		iter++;
	}
	
	return 1;
}
*/

/*double functionSet::objective(double *Z, double *Zstar, double *Q, double *SLOPE)
{

	error = 0;
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
		if(pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) > tetaC)
			error +=  	((-Z[j*M + i] + Zstar[j*M + i])/4 - D*(Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] -4 *Z[j*M + i])/(dx*dx) 
			+ K*(pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) - tetaC))
					*((-Z[j*M + i] + Zstar[j*M + i])/4 - D*(Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] -4 *Z[j*M + i])/(dx*dx) 
			+ K*(pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) - tetaC));
		else
			error +=  ((-Z[j*M + i] + Zstar[j*M + i])/4 - D*(Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] -4 *Z[j*M + i])/(dx*dx))
				*((Z[j*M + i] - Zstar[j*M + i])/4 - D*(Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] -4 *Z[j*M + i])/(dx*dx));
	
	return sqrt(error);
	
}
*/
/*int functionSet::linearModel(double *Z, double *Zstar,double *Q, double *SLOPE)
{

	double *XX, *YY;
	
	XX       = (double*)calloc(M*N*3 , sizeof(double)); //The design matrix for the linear model (XX[:,1] = diffusion; XX[:,2] = advection; XX[:,3] = 1) 
	YY       = (double*)calloc(M*N , sizeof(double));   //The experimental results in the linear model (time derivative) 

	
		
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		XX[0*M*N + (N-1)*M + i] = 0.0;
		XX[1*M*N + (N-1)*M + i] = 0.0;
		XX[2*M*N + (N-1)*M + i] = 0;
		YY[0*M*N + (N-1)*M + i] = 0;

		XX[0*M*N + j*M + M-1] = 0.0;
		XX[1*M*N + j*M + M-1] = 0.0;
		XX[2*M*N + j*M + M-1] = 0;
		YY[0*M*N + j*M + M-1] = 0;

		XX[0*M*N + j*M + i] = 	(
					(Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] -4 *Z[j*M + i])/(dx*dx) 
						+
					(Zstar[j*M + i + 1] + Zstar[j*M + i - 1] + Zstar[(j + 1)*M + i] + Zstar[(j - 1)*M + i] -4 *Zstar[j*M + i])/(dx*dx) 
					)/2;

		if(pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) > tetaC)
			XX[1*M*N + j*M + i] = pow(SLOPE[j*M + i], n)*pow(Q[j*M + i], m) - tetaC;
		else
			XX[1*M*N + j*M + i] = 0.0;

		XX[2*M*N + j*M + i] = 1.0;

		YY[0*M*N + j*M + i] = (-Z[j*M + i] + Zstar[j*M + i])/4.0;


	}

}

*/
