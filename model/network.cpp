
#include "header.h"


int functionSet::allocate(double *X,double *Y,double *Z,double *Q,double *dZstar, double *Zstar,double *PHI)
{
	int i,j;
	//X       = (double*)calloc(M*N , sizeof(double));
	param       = (double*)calloc(M*N , sizeof(double));
	CURVE  	= (double*)calloc(M*N , sizeof(double));
	R       = (double*)calloc(M*N , sizeof(double));
	D       = (double*)calloc(M*N , sizeof(double));
	alpha   = (double*)calloc(M*N , sizeof(double));
	Omega	= (double*)calloc(M*N , sizeof(double));
	K	= (double*)calloc(M*N , sizeof(double));
	
	DIR     = (int *) calloc(M*N, sizeof(int));
	INPUT   = (int *) calloc(M*N, sizeof(int));
	drainTo   = (int *) calloc(M*N, sizeof(int));
	drainFrom   = (int *) calloc(M*N, sizeof(int));
	//Q    = (double*)calloc(M*N , sizeof(double));

	//dZstar  = (double*)calloc(M*N , sizeof(double));
	//Zstar	= (double*)calloc(M*N , sizeof(double));
	//PHI   	= (double*)calloc(M*N , sizeof(double));
	a  = (double*)calloc(M*N , sizeof(double));
	b  = (double*)calloc(M*N , sizeof(double));
	c  = (double*)calloc(M*N , sizeof(double));
	d  = (double*)calloc(M*N , sizeof(double));

	down    = (double *) calloc(20*M, sizeof(double));
	up    = (double *) calloc(20*M, sizeof(double));
	left    = (double *) calloc(20*N, sizeof(double));
	right    = (double *) calloc(20*N, sizeof(double));
	scanTime    = (double *) calloc(20, sizeof(double));
	
	downBC    = (double *) calloc(M, sizeof(double));
	upBC    = (double *) calloc(M, sizeof(double));
	leftBC    = (double *) calloc(N, sizeof(double));
	rightBC    = (double *) calloc(N, sizeof(double));

	
	for(j=0; j<N; j++)
	for(i=0; i<M ; i++)
	{
		X[j*M + i] = 0; Y[j*M + i] = 0; Z[j*M + i] = 0; DIR[j*M + i] = 0; INPUT[j*M + i] = 0; drainTo[j*M + i] = 0;
		Q[j*M + i] = 0; dZstar[j*M + i] = 0; Zstar[j*M + i] = 0;PHI[j*M + i] = 0, CURVE[j*M + i] = 0;
	}
	
	return 0;
}


void functionSet::defineRain(string str1, int uni)//gets the normalized rainfall distribution and calculates the rainfall based ont the avergae rainfall given in the header file (aveRain)
{

	int i,j;
double number;
	input1.close();
	input1.open(str1.c_str());//
	if (uni ==0)
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
			R[j*M+i] = 1.0;
	}
	else
	{
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
		for(j=0; j<N; j++)	
			for(i=0; i<M; i++)
	    		{
			
				input1>>number;
				input1>>number;
				input1>>R[j*M+i]; 
				if(R[j*M+i] < 0)
					R[j*M+i] = 0;
				R[j*M+i] *= aveRain;        //for variable rainfall it is R[j*M+i] = *aveRain;
				//cout << X[j*M+i] << "    " <<Y[j*M+i]<< "    "<<Z[j*M+i]<<"    "<< endl;
	    		}
	


		input1.close();  
		}
	}
/*
	for(j=0; j<int(500.0/dy) ; j++)
	{
		for(i=0; i<int(200.0/dx); i++)
			R[j*M+i] = 0.0;
		
		for(i=int(200.0/dx); i<M; i++)
			R[j*M+i] = 0.0;
	}

	for(j=int(500.0/dy); j<int(700.0/dy); j++)
	{
		for(i=0; i<int(200.0/dx); i++)
			R[j*M+i] = 0.0;
		
		for(i=int(200.0/dx); i<M; i++)
			R[j*M+i] = 1.0;
	}

	for(j=int(700.0/dy); j<N; j++)
	{
		for(i=0; i<int(200/dx); i++)
			R[j*M+i] = 0.0;
		
		for(i=int(200/dx); i<M; i++)
			R[j*M+i] = 0.0;
	}
*/
}
	
void functionSet::defineD(int varD)
{
	int i,j;
	if(varD == 0)
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
			D[j*M+i] = Dave;//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave
	}
	else 
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
			D[j*M+i] = Dave*R[j*M+i]/aveRain;//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave	
	}

	
}


void functionSet::defineRandomD(int varD) // 0 or 1 is put directly in the function
{
	srand ( time(NULL) );
	int i,j;
	double random;
	std::random_device rd;
    	std::default_random_engine generator(rd());
	std::cauchy_distribution<double> cauchy(0,0.11757900431008157);
	if(varD == 0)
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
			D[j*M+i] = Dave;//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave
	}
	else 
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
		{
			random = 10;
			while (random>1)
				random = abs(cauchy(generator));
			
			D[j*M+i] = Dave*(1-random);//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave	
		}
	}	
	
}


void functionSet::defineK(int varK) //varK ==0:fixed K, varK ==1:random K
{
	int i,j;
	double random;
	std::random_device rd;
    	std::default_random_engine generator(rd());
	std::cauchy_distribution<double> cauchy(0,0.11757900431008157);
	if(varK == 0)
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
			K[j*M+i] = Kave;//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave
	}
	else 
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
		{
			random = 10;
			while (random>1)
				random = abs(cauchy(generator));
			K[j*M+i] = Kave*(1-random);//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave	

		}
	}	
}

void functionSet::defineOmega(int varOmega)
{
	int i,j;
	double random;
	std::random_device rd;
    	std::default_random_engine generator(rd());
	std::cauchy_distribution<double> cauchy(0,0.11757900431008157);
	if(varOmega == 0)
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
			Omega[j*M+i] = OmegaAve;//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave
	}
	else 
	{
		for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
		{
			random = 10;
			while (random>1)
				random = abs(cauchy(generator));
			Omega[j*M+i] = OmegaAve*random;//D[j*M+i] = Dave*R[j*M+i]/aveRain;D[j*M+i] = Dave
		}	
	}	
}

void functionSet::readLand(double *X, double *Y, double *Z, string str1)
{
	int i,j;
	double number;
	input1.close();
	input1.open(str1.c_str());//
	
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
		
//FOR RUNNING FROM THE BEGINING (ARCGIS DATA)	
	for(j=0; j<N; j++)	
		for(i=0; i<M; i++)
    		{
			
			input1>>X[j*M+i];
			input1>>Y[j*M+i];
			input1>>Z[j*M+i]; 
			//cout << X[j*M+i] << "    " <<Y[j*M+i]<< "    "<<Z[j*M+i]<<"    "<< endl;
    		}
	

//FOR RUNNING AFTER A CERTAIN TIME

			
	
	/*input1>>str1; 
	input1>>str1; 
	input1>>str1; 
	input1>>str1;
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	
    		{
			input1>>X[j*M+i];
			input1>>Y[j*M+i];
			input1>>Z[j*M+i]; 
			input1>>number;	
    		}	

		*/
	input1.close();  
	}

///FOR ANALYTICAL SOLUTION
	/*
	for(i=0; i<M; i++)
	for(j=0; j<N; j++)
	{
		X[j*M+i] = 2.0*i/(M - 1)-1.0; 
		Y[j*M+i] = 2.0*j/(N - 1)-1.0;
		Z[j*M+i] = 0.0;
		
	}
*/
	//cout << "readLand() is finished" << endl;
	
}



void functionSet::smoothZ(double *Z)
{
	int i,j;
	for(j=1; j<N-1; j++)
	for(i=1; i<M-1; i++)
		Z[j*M + i] = 0.25*Z[j*M + i] + 0.125*(Z[j*M + i + 1] + Z[(j + 1)*M + i] + Z[j*M + i - 1] + Z[(j - 1)*M + i])
					     + 0.0625*(Z[(j + 1)*M + i + 1] + Z[(j + 1)*M + i - 1] + Z[(j - 1)*M + i - 1] + Z[(j - 1)*M + i + 1]);
}




int functionSet::filling(double *Z, double *Zstar, double *Y)
{
	int i,j;
	double min, count = 0;

	for(i = 0; i < M; i++)
	{
		
		if(Z[0*M + i]<= Z[1*M + i]+epsilon)
		Z[0*M + i]= Z[1*M + i] +epsilon ;
		
		
		if(Z[(N - 1)*M + i]<= Z[(N - 2)*M + i]+epsilon )
			Z[(N - 1)*M + i]= Z[(N - 2)*M + i]+epsilon;
		
	}
	
	Z[0*M + M/2] = 0;	
	
	for(j = 0; j < N  ; j++)
	{

		
		if(Z[j*M + 0]<= Z[j*M + 1]+epsilon )
			Z[j*M + 0]= Z[j*M + 1] + epsilon;

		

		if(Z[j*M + M - 1]<= Z[j*M + M - 2]+epsilon )
			Z[j*M + M - 1]= Z[j*M + M - 2] + epsilon;
	}
	Z[0*M + 0] = (Z[1*M + 0] + Z[0*M + 1])/2;
	Z[0*M + M - 1] = (Z[1*M + M - 1] + Z[0*M + M - 2])/2;
	Z[(N - 1)*M + 0] = (Z[(N - 1)*M + 1] + Z[(N - 2)*M + 0])/2;
	Z[(N - 1)*M + M - 1] = (Z[(N - 2)*M +  M - 1] + Z[(N - 1)*M +  M - 2])/2;




	
	for(j=1; j<N-1; j++)
	{
		Zstar[j*M + 0] = Z[j*M + 0];
		Zstar[j*M + M - 1] = Z[j*M + M - 1];
	}

	for(i=1; i<M-1; i++)
	{
		Zstar[0*M + i] = Z[0*M + i];
		Zstar[(N-1)*M + i] = Z[(N-1)*M + i];
	}
	
	for(j=N-2; j>0; j--)
	for(i=1; i<M-1; i++)
	{
		Zstar[j*M + i] = 20000;
	}

	for(j=1; j<N-1; j++)
	for(i=1; i<M-1; i++)
	{
	
		min = Z[j*M + i + 1];
		if (Z[(j - 1)*M + i + 1] < min)
			min = Z[(j - 1)*M + i + 1];
		
		if (Z[(j - 1)*M + i] < min)
			min = Z[(j - 1)*M + i];

		if (Z[(j - 1)*M + i - 1] < min)
			min = Z[(j - 1)*M + i - 1];

		if (Z[(j - 0)*M + i - 1] < min)
			min = Z[(j - 0)*M + i - 1];

		if (Z[(j + 1)*M + i - 1] < min)
			min = Z[(j + 1)*M + i - 1];

		if (Z[(j + 1)*M + i - 0] < min)
			min = Z[(j + 1)*M + i - 0];

		if (Z[(j + 1)*M + i + 1] < min)
			min = Z[(j + 1)*M + i + 1];
		
		if (Z[j*M + i] > min + epsilon)
			Zstar[j*M + i] = Z[j*M + i];
		else if (Zstar[j*M + i] > min + epsilon )
		{
			Zstar[j*M + i] = min + epsilon ;
			count++;
		}
	}
	
	
	//cout <<"modified: "<< count << endl;
	for(j=0; j<N; j++)
	for(i=0; i<M; i++)
	Z[j*M + i] = Zstar[j*M + i];
	

	for(j=0; j<N; j++)
	for(i=0; i<M; i++)
		Zstar[j*M + i] = 0;
		
	return 0;
}
		
	
void functionSet::direction(double *Z, double *Zstar, double *Y)
{
	int i,j;
	double max;

	for(j=0; j<N; j++)
		for(i=0; i<M; i++)
			drainTo[j*M + i] = -1;
	for(j=1; j<N-1; j++)
	{
		DIR[j*M + 0] = 1;
		drainTo[j*M + 0] = j*M + 1;

		DIR[j*M + M - 1] = 16;
		drainTo[j*M + M - 1] = j*M + M - 2;
		
		/*DIR[j*M + 0] = -1;
		drainTo[j*M + 0] = -1;

		DIR[j*M + M - 1] = -1;
		drainTo[j*M + M - 1] = -1;*/
	}

	for(i=1; i<M-1; i++)
	{
		DIR[0*M + i] = 64;
		drainTo[0*M + i] = 1*M + i;

		DIR[(N-1)*M + i] = 4;
		drainTo[(N-1)*M + i] = (N-2)*M + i;


	
		/*DIR[0*M + i] = -1;
		drainTo[0*M + i] = -1;

		DIR[(N-1)*M + i] = -1;
		drainTo[(N-1)*M + i] = -1;*/

	}
	/*for(i = M/2-12; i <= M/2 + 12; i++)
	{
		DIR[0*M + i] = -1;
		drainTo[0*M + i] = -1;
	}
	i = 0;j = 0;
	DIR[i*M + i] = 7;
	drainTo[0*M + i] = (j + 1)*M + i + 1;

	i = M-1;j = 0;
	DIR[i*M + i] = 5;
	drainTo[0*M + i] = (j + 1)*M + i - 1;

	i = M-1;j = N-1;
	DIR[i*M + i] = 3;
	drainTo[0*M + i] = (j - 1)*M + i - 1;


	i = 0;j = N-1;
	DIR[i*M + i] = 1;
	drainTo[0*M + i] = (j - 1)*M + i + 1;
*/
	i = M/2;	
	DIR[0*M + i] = -1;
	drainTo[0*M + i] = -1;
	DIR[1*M + i] = 4;
	drainTo[1*M + i] = 0*M + i;

	double count = 10, count1 = 0;
	while (count>0)
	{
		int maxIndex;
		for(j=1; j<N-1; j++)
		for(i=1; i<M-1; i++)
		{

			slopeLocal[0] = (Z[j*M + i] - Z[j*M + i + 1]);
			slopeLocal[1] = (Z[j*M + i] - Z[(j - 1)*M + i + 1])/1.4142;
			slopeLocal[2] = (Z[j*M + i] - Z[(j - 1)*M + i]);
			slopeLocal[3] = (Z[j*M + i] - Z[(j - 1)*M + i - 1])/1.4142;
			slopeLocal[4] = (Z[j*M + i] - Z[j*M + i - 1]);
			slopeLocal[5] = (Z[j*M + i] - Z[(j + 1)*M + i - 1])/1.4142;
			slopeLocal[6] = (Z[j*M + i] - Z[(j + 1)*M + i]);
			slopeLocal[7] = (Z[j*M + i] - Z[(j + 1)*M + i + 1])/1.4142;	

			maxIndex = 0;      
			max = slopeLocal[0] ;
			int k;
			for (k = 1; k < 8; k++)
				if (slopeLocal[k] > max)
				{
					max = slopeLocal[k];
					maxIndex = k;
				}

			DIR[j*M + i] = -1;
			drainTo[j*M + i] = -1;

			if (maxIndex == 0 )
			{
				drainTo[j*M + i] = j*M + i + 1;
				DIR[j*M + i] = 1;
			}
			if (maxIndex == 1 )
			{
				drainTo[j*M + i] = (j - 1)*M + i + 1;
				DIR[j*M + i] = 2;
			}
			if (maxIndex == 2 )
			{
				drainTo[j*M + i] = (j - 1)*M + i;
				DIR[j*M + i] = 4;
			}
			if (maxIndex == 3 )
			{
				drainTo[j*M + i] = (j - 1)*M + i - 1;
				DIR[j*M + i] = 8;
			}
			if (maxIndex == 4 )
			{
				drainTo[j*M + i] = j*M + i - 1;
				DIR[j*M + i] = 16;
			}
			if (maxIndex == 5)
			{
				drainTo[j*M + i] = (j + 1)*M + i - 1;
				DIR[j*M + i] = 32;
			}
			if (maxIndex == 6)
			{
				drainTo[j*M + i] = (j + 1)*M + i;
				DIR[j*M + i] = 64;
			}
			if (maxIndex == 7)
			{
				drainTo[j*M + i] = (j + 1)*M + i + 1;
				DIR[j*M + i] = 128;
			}

		}
		count1 = 0;
		for(j=1; j<N-1; j++)
		for(i=1; i<M-1; i++)
		 	if(drainTo[drainTo[j*M + i]] == j*M + i || drainTo[drainTo[drainTo[j*M + i]]] == j*M + i || drainTo[drainTo[drainTo[drainTo[j*M + i]]]] == j*M + i || drainTo[j*M + i] == -1)
				count1++;
		
		if (count1 > 0)
			filling(Z, Zstar, Y);
		//cout << "PIT POIT: " << count1 << endl;	
		count = count1;
	}
	//cout << "direction() is finished" << endl;
}
		
	

void functionSet::drainageQ(double *Q)
{
	int i,j;
	for(j=0; j<N; j++)
	{
		Q[j*M + 0] = R[j*M + 0] ;
		INPUT[j*M + 0] = 0;

		Q[j*M + M - 1] = R[j*M + M - 1];
		INPUT[j*M + M - 1] = 0;
		
	}

	for(i=0; i<M; i++)
	{
		Q[0*M + i] = R[0*M + i] ;
		INPUT[0*M + i] = 0;

		Q[(N-1)*M + i] = R[(N-1)*M + i];
		INPUT[(N-1)*M + i] = 0;
	}
	
	for(j=N-2; j>0; j--)
	for(i=1; i<M-1; i++)
	{
		Q[j*M + i] = R[j*M + i];
		INPUT[j*M + i] = 0;	
	
		if (DIR[j*M + i + 1] == 16)
			INPUT[j*M + i] += 1;
	
		if (DIR[(j + 1)*M + i + 1] == 8)
			INPUT[j*M + i] += 1;
	
		if (DIR[(j + 1)*M + i] == 4)
			INPUT[j*M + i] += 1;

		if (DIR[(j + 1)*M + i - 1] == 2)
			INPUT[j*M + i] += 1;

		if (DIR[j*M + i - 1] == 1)
			INPUT[j*M + i] += 1;

		if (DIR[(j - 1)*M + i - 1] == 128)
			INPUT[j*M + i] += 1;

		if (DIR[(j - 1)*M + i] == 64)
			INPUT[j*M + i] += 1;

		if (DIR[(j - 1)*M + i + 1] == 32)
			INPUT[j*M + i] += 1;
		
	}


	
	int count = 10;
	int count1;

	while (count > 0)
	{
		count1=0;
		for(j=N-1; j>-1; j--)
		for(i=0; i<M; i++)
		{
			if(INPUT[j*M + i] == 0)
			{
				count1++;				
			 	Q[drainTo[j*M + i]] += Q[j*M + i];
				INPUT[drainTo[j*M + i]] -= 1;
				INPUT[j*M + i] = 1000;
			}
		}
	
		count = count1;
		
	}

	count1=0;
	for(j=N-2; j>0; j--)
	for(i=1; i<M-1; i++)
	{
		if(INPUT[j*M + i] != 1000)
		{
			count1++;				
			//cout << j*M + i<<"   input  "<< INPUT[j*M + i] << " Q  "<<Q[j*M + i]<< endl;	
		}
		
	}
	
	//cout << "PIT POIT2: " << count1 << endl;	
		
	for(j=0; j<N; j++)
	for(i=0; i<M; i++)
		Q[j*M + i] *=dx*dy;
	//cout << "drainageQ() is running" << endl;
	

}



void functionSet::slopeCalc(double *Z, double *SLOPE)
{
	double s[4];
	int i,j;
	
	
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
	{
		s[0] = (Z[j*M + i + 1] - Z[j*M + i - 1])/(2*dx);
		s[1] = (Z[(j + 1)*M + i] - Z[(j - 1)*M + i])/(2*dy);
		s[2] = (Z[(j - 1)*M + i + 1] - Z[(j + 1)*M + i - 1])/(2*sqrt(dx*dx + dy*dy));
		s[3] = (Z[(j + 1)*M + i + 1] - Z[(j - 1)*M + i - 1])/(2*sqrt(dx*dx + dy*dy));
		SLOPE[j*M + i] =  (sqrt(s[0]*s[0]+ s[1]*s[1]) + sqrt(s[2]*s[2]+ s[3]*s[3]))/2.0;
	}
	for(j=0; j<N; j++)
	{
		SLOPE[j*M + 0] =SLOPE[j*M + 1];SLOPE[j*M + M - 1] =SLOPE[j*M + M - 2];
	}

	for(i=0; i<M; i++)
	{
		SLOPE[(N-1)*M + i] =SLOPE[(N-2)*M + i];SLOPE[0*M + i] =SLOPE[1*M + i];
	}
	/*int i,j;
	for(j=0; j<N; j++)
		SLOPE[j*M + 0] =1e-15;SLOPE[j*M + M - 1] =1e-15;
	for(i=0; i<M; i++)
		SLOPE[(N-1)*M + i] =1e-15;SLOPE[0*M + i] =1e-15;
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
	{
		SLOPE[j*M + i] = (Z[j*M + i] - Z[drainTo[j*M + i]])/dx;

		if(DIR[j*M + i] == 2 || DIR[j*M + i] == 8 || DIR[j*M + i] == 32 || DIR[j*M + i] == 128)
			SLOPE[j*M + i] /= sqrt(2);
	}*/




}

void functionSet::curveCalc(double *Z, double *CURVE)
{
	int i,j;
	
	
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
	{
		
		CURVE[j*M + i] = (Z[j*M + i + 1] + Z[j*M + i - 1] + Z[(j + 1)*M + i] + Z[(j - 1)*M + i] - 4*Z[j*M + i])/(dx*dx);
	}
	for(j=0; j<N; j++)
	{
		CURVE[j*M + 0] =CURVE[j*M + 1];CURVE[j*M + M - 1] =CURVE[j*M + M - 2];
	}

	for(i=0; i<M; i++)
	{
		CURVE[(N-1)*M + i] =CURVE[(N-2)*M + i];CURVE[0*M + i] =CURVE[1*M + i];
	}
	
}



