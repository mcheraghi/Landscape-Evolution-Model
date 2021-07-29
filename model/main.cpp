#include "header.h"

/*
The output folder, the input DEM file and the rain contour has to be specified,
*/
void mainRUN()
{
	
	double *X,*Y,*Z,*Q, *Qx, *Qy, *AREA,*SLOPE,*CURVE,*L,*LL,*downL,*Zstar,*dZstar,*PHI,*order; // the parameter dzstar in mainRUN is replaced by L
	int RunOpt, i, j;
	
	X       = (double*)calloc(M*N , sizeof(double));
	Y       = (double*)calloc(M*N , sizeof(double));
	Z       = (double*)calloc(M*N , sizeof(double));
	Q       = (double*)calloc(M*N , sizeof(double));
	Qx       = (double*)calloc(M*N , sizeof(double));
	Qy       = (double*)calloc(M*N , sizeof(double));
	SLOPE       = (double*)calloc(M*N , sizeof(double));
	CURVE       = (double*)calloc(M*N , sizeof(double));
	Zstar       = (double*)calloc(M*N , sizeof(double));
	dZstar       = (double*)calloc(M*N , sizeof(double));
	PHI       = (double*)calloc(M*N , sizeof(double));
	L       = (double*)calloc(M*N , sizeof(double));
	AREA       = (double*)calloc(M*N , sizeof(double));
	LL       = (double*)calloc(M*N , sizeof(double));
	downL       = (double*)calloc(M*N , sizeof(double));
	order       = (double*)calloc(M*N , sizeof(double));
	
	interface III;
	functionSet MCH;

	MCH.allocate(X,Y,Z,Q,dZstar, Zstar,PHI);
	  



	MCH.Dave =17571.0;	
	MCH.K = 0.18997;			
	MCH.m = 0.5;
	MCH.n = 1.0;
	MCH.tetaC = 13.0765;
	MCH.E = 0;
	MCH.rndCoef =0;
	char *path1=NULL;
    	size_t size;
    	string path=getcwd(path1,size);
    	cout<<"\n current Path"<<path;
  	path = path.substr(0,path.size()-6);     // take the "/model" part
	string addressRUN = path+"/results";

	//boost::filesystem::path full_path( boost::filesystem::current_path() );
	//full_path-="/model";
	
	//string addressRUN = "/home/mohsen/Man2/DATA/RUNS/run1/results"; //the output file
	MCH.varD = 0;	//0: uniform, 1: variable
		
	double initialTime = 0.25;//NOT 0.5 and 1.0! 
	double finalTime = 16.0;//NOT 0.5 and 1.0! 
	string initialTimeZ, finalTimeZ;

	if (dx == 8)
	{
		if(initialTime == 0.25)
			initialTimeZ = path+"/input/00h15min_8mm.dat"; //the input file
		else if(initialTime == 0.5)
			initialTimeZ = path+"/input/00h30min_8mm.dat"; //the input file
		else if(initialTime == 1.0)
			initialTimeZ = path+"/input/01h00min_8mm.dat"; //the input file
		else if(initialTime == 2.0)
			initialTimeZ = path+"/input/02h00min_8mm.dat"; //the input file
		else if(initialTime == 4.0)
			initialTimeZ = path+"/input/04h00min_8mm.dat"; //the input file
		else if(initialTime == 8.0)
			initialTimeZ = path+"/input/08h00min_8mm.dat"; //the input file
	
	
		if(finalTime == 0.25)
			finalTimeZ = path+"/input/00h15min_8mm.dat"; //the input file
		else if(finalTime == 2.0)
			finalTimeZ = path+"/input/02h00min_8mm.dat"; //the input file
		else if(finalTime == 4.0)
			 finalTimeZ = path+"/input/04h00min_8mm.dat"; //the input file
		else if(finalTime == 8.0)
			finalTimeZ = path+"/input/08h00min_8mm.dat"; //the input file
		else if(finalTime == 16.0)
			finalTimeZ = path+"/input/16h00min_8mm.dat"; //the input file

		MCH.defineRain(path+"/input/Rain8mm.dat",1); /// const = 0:uniform rain
	}
	
	else
	{

		
	if(initialTime == 0.25)
		initialTimeZ = "/home/mohsen/LandScapeCode/input/00h15min_4mm.dat"; //the input file
	if(initialTime == 0.5)
		initialTimeZ = "/home/mohsen/LandScapeCode/input/00h30min_4mm.dat"; //the input file
	if(initialTime == 1.0)
		initialTimeZ = "/home/mohsen/LandScapeCode/input/01h00min_4mm.dat"; //the input file
	else if(initialTime == 2.0)
		initialTimeZ = "/home/mohsen/LandScapeCode/input/02h00min_4mm.dat"; //the input file
	else if(initialTime == 4.0)
		initialTimeZ = "/home/mohsen/LandScapeCode/input/04h00min_4mm.dat"; //the input file
	else if(initialTime == 8.0)
		initialTimeZ = "/home/mohsen/LandScapeCode/input/08h00min_4mm.dat"; //the input file

	
	if(finalTime == 0.25)
		finalTimeZ = "/home/mohsen/LandScapeCode/input/00h15min_4mm.dat"; //the input file
	else if(finalTime == 2.0)
		finalTimeZ = "/home/mohsen/LandScapeCode/input/02h00min_4mm.dat"; //the input file
	else if(finalTime == 4.0)
		 finalTimeZ = "/home/mohsen/LandScapeCode/input/04h00min_4mm.dat"; //the input file
	else if(finalTime == 8.0)
		finalTimeZ = "/home/mohsen/LandScapeCode/input/08h00min_4mm.dat"; //the input file
	else if(finalTime == 16.0)
		finalTimeZ = "/home/mohsen/LandScapeCode/input/16h00min_4mm.dat"; //the input file

	MCH.defineRain("/home/mohsen/LandScapeCode/input/Rain4mm.dat",1); ///
	}	
	string filename[] = {path+"/input/00h15min_8mm.dat", path+"/input/00h30min_8mm.dat", path+"/input/01h00min_8mm.dat", path+"/input/02h00min_8mm.dat", path+"/input/04h00min_8mm.dat", path+"/input/08h00min_8mm.dat", path+"/input/16h00min_8mm.dat"};
	
	if (dx == 4) //attention that the file names have the same size
	{
		string filename[] = {"/home/mohsen/LandScapeCode/input/00h15min_4mm.dat", "/home/mohsen/LandScapeCode/input/00h30min_4mm.dat", "/home/mohsen/LandScapeCode/input/01h00min_4mm.dat","/home/mohsen/LandScapeCode/input/02h00min_4mm.dat", "/home/mohsen/LandScapeCode/input/04h00min_4mm.dat", "/home/mohsen/LandScapeCode/input/08h00min_4mm.dat", "/home/mohsen/LandScapeCode/input/16h00min_4mm.dat"};
		cout <<filename[1] << endl;
 
	}
	MCH.temporalBCset(X, Y, Z,filename);// captures the BCs from the scans
	MCH.RUN(X, Y, Z, Zstar, Q, Qx, Qy, SLOPE, CURVE, PHI, dZstar, AREA, LL,downL, addressRUN, initialTime, finalTime, initialTimeZ,  finalTimeZ);//0:constant D, 1: variable D [refer to defineD(int varD) function]
	

	//MCH.mainRiver(X, Y, Z, Zstar, Q, SLOPE, CURVE, PHI, dZstar, addressRUN);
	
}


void mainOPT()
{
	time_t tstart, tend; 
 	tstart = time(0);
	double *X,*Y,*Z,*Q,*SLOPE,*CURVE,*dZstar,*Zstar,*PHI;
	
	
	X       = (double*)calloc(M*N , sizeof(double));
	Y       = (double*)calloc(M*N , sizeof(double));
	Z       = (double*)calloc(M*N , sizeof(double));
	Q       = (double*)calloc(M*N , sizeof(double));
	SLOPE       = (double*)calloc(M*N , sizeof(double));
	CURVE       = (double*)calloc(M*N , sizeof(double));
	dZstar       = (double*)calloc(M*N , sizeof(double));
	Zstar       = (double*)calloc(M*N , sizeof(double));
	PHI       = (double*)calloc(M*N , sizeof(double));
	
	interface III;
	functionSet MCH;
	double T1 = 0.25;
	double T2 = 8;
	
	char *path1=NULL;
    	size_t size;
    	string path=getcwd(path1,size);
    	cout<<"\n current Path"<<path;
  	path = path.substr(0,path.size()-6);     // take the "/model" part
	string addressOPT = path+"/results";


	MCH.allocate(X,Y,Z,Q,dZstar, Zstar,PHI);
	//string addressOPT = "/home/mohsen/LandScapeCode/2D/Results/exp8/OPTIMIZATION/interpolatedBC/confidenseInterval/uniD_opAt8h.dat";	 //the output file for optimization
	MCH.varD = 0;
	string initialTime = path+"/input/00h15min_8mm.dat"; //the input file
	string finalTime = path+"/input/08h00min_8mm.dat"; //the input file
	MCH.defineRain(path+"/input/Rain8mm.dat",0); ///
	
		string filename[] = {path+"/input/00h15min_8mm.dat", path+"/input/00h30min_8mm.dat", path+"/input/01h00min_8mm.dat", path+"/input/02h00min_8mm.dat", path+"/input/04h00min_8mm.dat", path+"/input/08h00min_8mm.dat", path+"/input/16h00min_8mm.dat"};
	
	if (dx == 4) //attention that the file names have the same size
	{
		string filename[] = {"/home/mohsen/LandScapeCode/input/00h15min_4mm.dat", "/home/mohsen/LandScapeCode/input/00h30min_4mm.dat", "/home/mohsen/LandScapeCode/input/01h00min_4mm.dat","/home/mohsen/LandScapeCode/input/02h00min_4mm.dat", "/home/mohsen/LandScapeCode/input/04h00min_4mm.dat", "/home/mohsen/LandScapeCode/input/08h00min_4mm.dat", "/home/mohsen/LandScapeCode/input/16h00min_4mm.dat"};
		cout <<filename[1] << endl;
 
	}
	MCH.temporalBCset(X, Y, Z,filename);// captures the BCs from the scans
	MCH.calibratePSO(X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, addressOPT, initialTime,finalTime,filename);


	double q[3] = {17571,0.184997,13.0765};
	double percentage[3] = {0.1,0.1,0.1};
	//MCH.confidenceInterval1(q, percentage, X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1,T2,initialTime, finalTime);
	//MCH.confidenceInterval(q, percentage, X, Y, Z, Zstar, Q, SLOPE, PHI, dZstar, T1,T2,initialTime, finalTime);

	/*X[0*2 + 0] = 25.1;
	X[1*2 + 0] = 2.2;
	X[2*2 + 0] = 4.1;
	X[0*2 + 1] = 3;
	X[1*2 + 1] = 5;
	X[2*2 + 1] = 6;

	MCH.transpose(X, Y, 2, 3);
	
	MCH.multiple(Y, X, Z, 3, 2, 2, 3);
	for(int i  = 0; i < 9; i++)
	{
		cout << Z[i]<<endl;
	}*/

	tend = time(0); 
	cout << "It took "<< difftime(tend, tstart)/3600 <<" hours"<< endl;
	
}

void mainAnalysis()
{
	
	double *X,*Y,*Z,*Zmodel,*Q, *Qx, *Qy, *AREA,*SLOPE,*CURVE,*downL,*L,*LL,*Zstar,*PHI, *power, *SLOPE_CURVE, *order; // the parameter dzstar in mainRUN is replaced by L
	int RunOpt, i, j;
	
	X       = (double*)calloc(M*N , sizeof(double));
	Y       = (double*)calloc(M*N , sizeof(double));
	Z       = (double*)calloc(M*N , sizeof(double));
	Zmodel       = (double*)calloc(M*N , sizeof(double));
	Q       = (double*)calloc(M*N , sizeof(double));
	Qx       = (double*)calloc(M*N , sizeof(double));
	Qy       = (double*)calloc(M*N , sizeof(double));
	SLOPE       = (double*)calloc(M*N , sizeof(double));
	CURVE       = (double*)calloc(M*N , sizeof(double));
	SLOPE_CURVE = (double*)calloc(M*N , sizeof(double));
	Zstar       = (double*)calloc(M*N , sizeof(double));
	PHI       = (double*)calloc(M*N , sizeof(double));
	L       = (double*)calloc(M*N , sizeof(double));
	AREA       = (double*)calloc(M*N , sizeof(double));
	LL       = (double*)calloc(M*N , sizeof(double));
	downL       = (double*)calloc(M*N , sizeof(double));
	power       = (double*)calloc(M*N , sizeof(double));
	order      = (double*)calloc(M*N , sizeof(double));
	
	interface      III;
	functionSet MCH;

	MCH.allocate(X,Y,Z,Q,L, Zstar,PHI);
	MCH.epsilon = 0.05;
	string r = "0.5";


	char *path1=NULL;
    	size_t size;
    	string path=getcwd(path1,size);
    	cout<<"\n current Path"<<path;
  	path = path.substr(0,path.size()-6);     // take the "/model" part
	string  folder = path+"/results";


	string fileName[] = {   path+"/input/00h15min_8mm.dat", 
				path+"/input/00h30min_8mm.dat", 
				path+"/input/01h00min_8mm.dat", 
				path+"/input/02h00min_8mm.dat", 
				path+"/input/04h00min_8mm.dat", 
				path+"/input/08h00min_8mm.dat", 
				path+"/input/16h00min_8mm.dat"};
				
	MCH.defineRain(path+"/input/Rain8mm.dat",0); ///
	string modelName[] = {"/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ00..dat", "/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ01..dat", "/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ02..dat", "/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ03..dat", "/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ04..dat", "/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ05..dat", "/home/mohsen/LandScapeCode/2D/Results/exp8/RUNS/87_8mm/XYZ06..dat"};



	string outName[] = {"0.25h.csv", "0.5h.csv", "1h.csv", "2h.csv", "4h.csv", "8h.csv", "16h.csv"}; 

	int ii, jj;
	for (i = 0; i < 7; i++)
	{
	
		string input = fileName[i]; //the input file
		string makefolder = "mkdir " + folder;
		system(makefolder.c_str());
		MCH.readLand(X, Y, Z,input);
		//MCH.filternn(Z,Zstar, 30);
		//MCH.replaceError(Z,Zstar); //Zstar to Z
		MCH.defineRain("/home/mohsen/LandScapeCode/input/Rain8mm.dat",0);
		MCH.temporalBCset(X, Y, Zstar,fileName);
		MCH.filling(Z, PHI, Y);
		MCH.direction(Z, PHI, Y);
		
		MCH.drainageQ(Q);
		MCH.drainageAREAandLL(AREA, LL,order);
		MCH.L_down(downL);
		MCH.slopeCalc(Z,SLOPE);
		MCH.curveCalc(Z,CURVE);
		//MCH.Qcomponents(X, Y, Q, Qx,Qy);
		MCH.writeDATA(X, Y,Z,Z,SLOPE,Q,AREA,LL, CURVE,downL, folder, i);
		MCH.writeXYZ(X, Y,Z, folder, i);

}
		/*
		int count = 0;
		for(ii = 0; ii<M; ii++)
		for(jj = 0; jj<N; jj++)
		{
			if(SLOPE[jj*M + ii]>0)
			{

				Y[count] =SLOPE[jj*M + ii];
				X[count] = AREA[jj*M + ii];
				count++;
			}
			
			if(SLOPE[jj*M + ii]>0 && CURVE[jj*M + ii]>0)
			{
				Zstar[jj*M + ii] = SLOPE[jj*M + ii]/CURVE[jj*M + ii];
				Y[count] =SLOPE[jj*M + ii]/CURVE[jj*M + ii];
				X[count] = AREA[jj*M + ii];
				count++;	
			}
		}
		//MCH.Kmean(Q,SLOPE, PHI,20, folder, "/16h_Kmeans_Q_SLOPE1.ods");//
		//MCH.Kmean(AREA,SLOPE, PHI,20, folder, "/16h_Kmeans_AREA_SLOPE1.ods");//
		//MCH.Kmean(Q,Zstar, PHI,20, folder, "/16h_Kmeans_Q_SLOPEstar1.ods");//
		//MCH.Kmean(X,Y, PHI,M*N,20, folder, "/Kmeans_AREA_SLOPE"+outName[i]);//
		
		ofstream out;
		out.close();
		string folderName = folder + "/AREA_SLOPE"+outName[i];
		remove(folderName.c_str());
		out.open(folderName.c_str(),ios::app);
		out << "AREA" <<",   " << "SLOPE" << endl;
		for(ii = 0; ii<M; ii++)
		for(jj = 0; jj<N; jj++)
			out << AREA[jj*M + ii] <<",   " << SLOPE[jj*M + ii] << endl;
		*/


		//MCH.curveCalc(Z,CURVE);
		//MCH.filter(CURVE,Zmodel,1);
		//MCH.filter(CURVE,Q,2);
		//MCH.filter(CURVE,Qx,4);
		//MCH.filter(CURVE,Qy,8);
		//MCH.filter(CURVE,SLOPE,16);
		//MCH.writeDATA(X, Y, Z, CURVE,Zmodel,Q,Qx,Qy, SLOPE, folder, i);
	
	

/*	//This calcualtes the data based on the morphology. it was used for the python figures in Mans2:LEM at the flume scales
	for (i = 0; i < 7; i++)
	{
		string input = fileName[i]; //the input file
		string inputModel = modelName[i]; //the input file
		string makefolder = "mkdir " + folder;
		system(makefolder.c_str());
	
	
		MCH.readLand(X, Y, Z,input);
		MCH.readLand(X, Y, Zmodel,inputModel);
	
	
		MCH.defineRain("/home/mohsen/LandScapeCode/input/Rain8mm.dat"); ///
	
		MCH.temporalBCset(X, Y, Zstar);// captures the BCs from the scans(don't put Z in this function)
	
		//MCH.boundarySet(Z, Y);
	
		//MCH.filterNetwork(Zmodel, Z, Y, Q);
	
		//MCH.filter(Z,Zstar,50);
		
		MCH.filling(Z, PHI, Y);
		MCH.direction(Z, PHI, Y);
		MCH.drainageQ(Q);
		MCH.slopeCalc(Z,SLOPE);
		MCH.filter(SLOPE,SLOPE,10);
		MCH.curveCalc(Z,CURVE);
		MCH.filter(CURVE,CURVE,1);	
	
		MCH.Qcomponents(X, Y, Q, Qx,Qy);
		//MCH.filter(Q,Q,1);
		//MCH.filterNetwork(Q,Q,3,1e6);
		//MCH.filter(Qx,Qx,4);
		//MCH.filter(Qy,Qy,4);
		MCH.writeDATA(X, Y, Z, SLOPE,Q, Qx, Qy, CURVE, folder, i);
	}
	
	MCH.drainageAREAandLL(AREA, LL);
	MCH.slopeCalc(Z,SLOPE);
	MCH.filter(SLOPE,SLOPE,10);
	MCH.curveCalc(Z,CURVE);
	MCH.filter(CURVE,CURVE,1);	
	MCH.writeDATA(X, Y, Z, SLOPE,Q,CURVE, folder, 8);
	MCH.filter(Q,Q,8);
	MCH.writeDATA(X, Y, Z, SLOPE,Q,CURVE, folder, 16);
	MCH.filter(Q,Q,16);
	MCH.writeDATA(X, Y, Z, SLOPE,Q,CURVE, folder, 32);*/

	//MCH.calcDrainFrom(Q);
	//MCH.UpLenght(L);//is based on the drainage area
	
	//MCH.mainRiver(X, Y, Z, L, Q, 2.0*aveRain*dx*dy, SLOPE, Q, CURVE, LL, PHI, folder, "/maxQriver.dat");
	
	//MCH.calcDrainFrom(AREA);
	//MCH.UpLenght(LL);//is based on the drainage area
	//MCH.mainRiver(X, Y, Z, L, AREA, 5.0*dy*dx, SLOPE, Q, CURVE, LL, PHI, folder,"/logestRiver.dat");
	

	//MCH.statistics(X, Y, Z, Q, AREA, SLOPE, CURVE, PHI, L, LL, 0.05, folder);//
	
	//MCH.Kmean(Q, SLOPE, PHI,100, folder, "/Kmean_K100_QS.txt");//
	//MCH.normalize(AREA);
	//MCH.normalize(SLOPE);
	//MCH.Ktest(AREA, SLOPE, PHI,20 , folder, "/KTest.txt");
	//MCH.Kmean(AREA, SLOPE, PHI, 50 , folder, "/Kmean_K50_AS.txt");
		
	
	
	/*for(j = 0; j<N; j++)
	for(i = 0; i<N; i++)
	{
		if(CURVE[j*M + i] !=0)
		SLOPE_CURVE[j*M + i] = SLOPE[j*M + i]/CURVE[j*M + i];
		else
		SLOPE_CURVE[j*M + i] = 0;
		

	}
	MCH.Kmean(Q, SLOPE_CURVE, PHI,50, folder, "/Kmean_K50_SLOPE_CURVEQ.txt");*/

	double min, max;
	MCH.findMinMax(AREA, &min,  &max);
	cout << min <<"    "<<max<<endl;

	MCH.findMinMax(Q, &min,  &max);
	cout << min <<"    "<<max<<endl;
}

int main()
{
	int mainNumb;
	
	cout << " HELLO Mohsen, Enter the main number (RUN:1, OPT:2 and ANALYSIS: 3)" << endl;
	cin >> mainNumb;
	cout << " Note: The resolution is: " <<dx<<" mm"<< endl;
	cout << " Je vais checker" << endl;
	if(mainNumb == 1)
	{
		cout << " OK, On y va :) !!!" << endl;
		mainRUN();
	}
	else if(mainNumb == 2)
	{
		cout << " OK, On y va :) !!!" << endl;
		mainOPT();
	}
	else if(mainNumb == 3)
	{
		cout << " OK, On y va :) !!!" << endl;
		mainAnalysis();
	}
	else
	{
		cout << "Desolé, Le shifre doit etre 1, 2, ou 3!!!!" << endl;
	}
	return 0;
}












/*ofstream output1;
	remove("/home/mohsen/LandScapeCode/input/dZdT15%80mmh_4mmGrd_12h.txt");
	output1.open("/home/mohsen/LandScapeCode/input/dZdT15%80mmh_4mmGrd_12h.txt",ios::app);

	MCH.readLand(X, Y, Z,"/home/mohsen/LandScapeCode/input/XYZ15%80mmh_4mmGrd_12h.txt");
	MCH.boundarySet(Z);
	MCH.filling(Z, Zstar);
	MCH.direction(Z, Zstar);
	
	MCH.readLand(X, Y, Q, "/home/mohsen/LandScapeCode/input/XYZ_15%80mmh_4mmGrd_20h.txt");
	MCH.boundarySet(Q);
	MCH.filling(Q, PHI);
	MCH.direction(Q, PHI);
	output1<<" X,    Y,    Z,    dzdt" <<endl;
	for(int i = 0; i<M; i++)
	for(int j = 0; j<N; j++)
	output1 << X[j*M + i] << ",    " << Y[j*M + i] <<",   " << Z[j*M + i] <<",   " << (Q[j*M + i] - Z[j*M + i])/8.0<< "   " << endl;

	*/





//for analztical solution, changes:dx dz delta,M,N, E = 1, K = 0, in PHI function if ( K==0 PHI = E),readinput function, 

	/*double ERROR = 5;
	int num = 0;
	double t = 0;
	double step = Tmax/10.0, captureTime = t + step;
	MCH.writeDATA(X, Y, Z, Q, num);
	while(t < Tmax)
	{
		
		MCH.dtCalc(Q, SLOPE);
		MCH.RungeKutta(Z, Zstar, Q, SLOPE, dZstar, PHI);
		MCH.crankNicholson(Z, Zstar);
		MCH.replaceError(Z, Zstar);
		MCH.boundarySet(Z);
		
		t+=MCH.dt;
		ERROR = MCH.error;
		cout <<"CHANGE (MM) =" << ERROR << "    " << "t =" << t << endl;
		if(t>captureTime)
		{MCH.writeDATA(X, Y, Z, Q, ++num);captureTime+=step;}

	}

	MCH.analytic(X, Y, Z);*/

