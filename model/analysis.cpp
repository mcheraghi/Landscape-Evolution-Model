#include "header.h"

void functionSet::drainageAREAandLL(double *AREA, double *LL, double *order) //the drainage area without considering the rainfall
{
	int i,j;
	//cout << "drainageQ() is running" << endl;
	for(j=0; j<N; j++)
	{
		AREA[j*M + 0] = 1 ;
		LL[j*M + 0] = dx;
		INPUT[j*M + 0] = 0;

		AREA[j*M + M - 1] = 1;
		LL[j*M + M - 1] = dx;
		INPUT[j*M + M - 1] = 0;
	}

	for(i=0; i<M; i++)
	{
		AREA[0*M + i] = 1 ;
		LL[0*M + i] = dy ;
		INPUT[0*M + i] = 0;

		AREA[(N-1)*M + i] = 1;
		LL[(N-1)*M + i] = dy;
		INPUT[(N-1)*M + i] = 0;
	}
	
	for(j=N-2; j>0; j--)
	for(i=1; i<M-1; i++)
	{
		AREA[j*M + i] = 1; // the initial area of each single cell
		INPUT[j*M + i] = 0;

		if(drainTo[j*M + i] == 	j*M + i + 1 || drainTo[j*M + i] == j*M + i - 1)// the initial lenght of each single cell
			LL[j*M + i] = dx;
		else if(drainTo[j*M + i] == (j - 1)*M + i || drainTo[j*M + i] == (j - 1)*M + i)
			LL[j*M + i] = dy;
		else
			LL[j*M + i] = sqrt(dx*dx + dy*dy);
	
	
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
	int mm, nn;
	int ord = 0;
	while (count > 0)
	{
		count1=0;
		
		for(j=N-1; j>-1; j--)
		for(i=0; i<M; i++)
		{

			if(INPUT[j*M + i] == 0)
			{
				INPUT[j*M + i] = ord;
				count1++;				
			 	AREA[drainTo[j*M + i]] += AREA[j*M + i];
				if(INPUT[drainTo[j*M + i]] == 1)
					LL[drainTo[j*M + i]] += LL[j*M + i];
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
		}
		
	}
	cout << "PIT POIT2: " << count1 << endl;	
		
	for(j=0; j<N; j++)
	for(i=0; i<M; i++)
		AREA[j*M + i] *=dx*dy;	

}

int functionSet::Qcomponents(double *X, double *Y, double *Q, double *Qx, double *Qy)
{
	double d;
	int i, j;
	for(j = 1; j < N - 1 ; j++)
	for(i = 1; i < M - 1; i++)	
	{
		d = sqrt((X[drainTo[j*M + i]] -  X[j*M + i])*(X[drainTo[j*M + i]] -  X[j*M + i]) + (Y[drainTo[j*M + i]] -  Y[j*M + i])*(Y[drainTo[j*M + i]] -  Y[j*M + i]));
		Qx[j*M + i] = Q[j*M + i]*(X[drainTo[j*M + i]] -  X[j*M + i])/d;
		Qy[j*M + i] = Q[j*M + i]*(Y[drainTo[j*M + i]] -  Y[j*M + i])/d;
	}
	
	
	for(j = 1; j < N - 1 ; j++)
	{
		i = 0;
		Qx[j*M + i] = Q[j*M + i];
		Qy[j*M + i] = 0;

		i = M - 1;
		Qx[j*M + i] = -Q[j*M + i];
		Qy[j*M + i] = 0;

	}
	

	for(i = 1; i < M - 1 ; i++)
	{
		j = 0;
		Qx[j*M + i] = 0;
		Qy[j*M + i] = Q[j*M + i];

		j = N - 1;
		Qx[j*M + i] = 0;
		Qy[j*M + i] = -Q[j*M + i];
	}

	
	i = M/2; j = 0;
	Qx[j*M + i] = 0;
	Qy[j*M + i] = -Q[j*M + i];


	
	Qx[0*M + 0] = (Qx[1*M + 0] + Qx[0*M + 1])/2;
	Qx[0*M + M - 1] = (Qx[1*M + M - 1] + Qx[0*M + M - 2])/2;
	Qx[(N - 1)*M + 0] = (Qx[(N - 1)*M + 1] + Qx[(N - 2)*M + 0])/2;
	Qx[(N - 1)*M + M - 1] = (Qx[(N - 2)*M +  M - 1] + Qx[(N - 1)*M +  M - 2])/2;

	Qy[0*M + 0] = (Qy[1*M + 0] + Qy[0*M + 1])/2;
	Qy[0*M + M - 1] = (Qy[1*M + M - 1] + Qy[0*M + M - 2])/2;
	Qy[(N - 1)*M + 0] = (Qy[(N - 1)*M + 1] + Qy[(N - 2)*M + 0])/2;
	Qy[(N - 1)*M + M - 1] = (Qy[(N - 2)*M +  M - 1] + Qy[(N - 1)*M +  M - 2])/2;
	
	return 1;
}

int functionSet::calcDrainFrom(double *BV)// instead of Q, the other parameters can be used (For example to find the maximum lenghts, the LL parameter can be used.
{
	int i,j;
// The grid should be filled and the directions be set (filling(), direction())before using this function 
	double max = 0;
	int maxPlace = 0;
	int mm, nn;
	for(j=0; j<N; j++)
	{
		drainFrom[j*M + 0] = j*M + 0;
		drainFrom[j*M + M - 1] = j*M + M - 1;
	}

	for(i=0; i<M; i++)
	{
		drainFrom[0*M + i] = 0*M + i ;
		drainFrom[(N-1)*M + i] = (N-1)*M + i;	
	}
	for(j = 1; j<N-1; j++)
	for(i = 1; i<M-1; i++)
	{
		max = 0;
		for(mm = -1; mm<2 ; mm++)
		for(nn = -1; nn<2 ; nn++)
		{
			if(drainTo[(j + nn)*M + i + mm] == j*M + i)
			{
				if(BV[(j + nn)*M + i + mm] > max)
				{
					max = BV[(j + nn)*M + i + mm];
					maxPlace = (j + nn)*M + i + mm;
				}
			}
		}
		drainFrom[j*M + i] = maxPlace ;
		
	}

}


int functionSet::mainRiver(double *X, double *Y, double *Z, double *L, double *BV, double BVmin,  double *SLOPE, double *Q, double *CURVE,double *LL, double *PHI, string folder, string fileName) //the main path is determined base on the based value (BV) and its minimum value (BVmin)
{
	int i,j;
	double max = 0;
	int maxPlace;
	int mm, nn;

	prepareFile(folder, fileName); //in prepares the output 1 for writing
	max = 0;
	for(i = 1; i<M-1; i++)
		if(BV[2*M + i] > max)
		{
			max = BV[2*M + i];
			maxPlace = 2*M + i;
		}
	output1 << "X" << "  " << "Y" << "   " << "Z" << "  " <<"L"<<"   "<<"LL"<<"   "<< "BaseValue" << "	" << "SLOPE" << "	"<< "DISCHARGE" <<"  "<<"CURVE" <<"  " "filteredCurve1" <<"  "<< "filteredCurve10" <<"  "<< "filteredSlope1" <<"  "<< "filteredSlope10" << endl;
		
	filter33(CURVE, a, 1); // a, b, c and d are used in the numerical process
	filter33(CURVE, b, 10);
	filter33(SLOPE, c, 1);
	filter33(SLOPE, d, 10);
		
	while(BV[maxPlace]>= BVmin)
	{
		output1 << X[maxPlace] << "	" << Y[maxPlace] << "	" << Z[maxPlace] <<"    "<<L[maxPlace] << "	"<<LL[maxPlace] << "	" << BV[maxPlace] << "	" << SLOPE[maxPlace] <<"	"<<Q[maxPlace] << "	"<<CURVE[maxPlace] <<"	 "<<a[maxPlace] <<"	"<<b[maxPlace] <<"	"<<c[maxPlace] <<"	"<<d[maxPlace] << endl;
		//cout << "maxPlace" <<endl;
		
			
		maxPlace = drainFrom[maxPlace];//drainFrom defines the input point with the maximum dischaerge
	}
			
	
	return 0;
}






void functionSet::UpLenght(double *L) // the length of the points are calculated based on the flow rate
{
	int i,j;
	for(j=0; j<N; j++)
	{
		L[j*M + 0] = dx ;
		INPUT[j*M + 0] = 0;

		L[j*M + M - 1] = dx;
		INPUT[j*M + M - 1] = 0;
		
	}

	for(i=0; i<M; i++)
	{
		L[0*M + i] = dy ;
		INPUT[0*M + i] = 0;

		L[(N-1)*M + i] = dy;
		INPUT[(N-1)*M + i] = 0;
	}
	
	
	for(j=N-2; j>0; j--)
	for(i=1; i<M-1; i++)
	{	
		if(drainTo[j*M + i] == 	j*M + i + 1 || drainTo[j*M + i] == j*M + i - 1)
			L[j*M + i] = dx;
		else if(drainTo[j*M + i] == (j - 1)*M + i || drainTo[j*M + i] == (j - 1)*M + i)
			L[j*M + i] = dy;
		else
			L[j*M + i] = sqrt(dx*dx + dy*dy);
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
		for(j=N-2; j>0; j--)
		for(i=1; i<M-1; i++)
		{
			if(INPUT[drainFrom[j*M + i]] == 0)
			{
				count1++;				
			 	L[j*M + i] += L[drainFrom[j*M + i]];
				INPUT[j*M + i] = 0;
				INPUT[drainFrom[j*M + i]] = 1000;
			}
		}
	
		count = count1;
		//cout<<count1<<endl;
		
	}
	
}

int functionSet::statistics(double *X, double *Y, double *Z, double *Q, double *AREA, double *SLOPE, double *CURVE, double *PHI, double *L, double *LL, double percentage, string folder)
//percentage(0-1) defines the ratio of points to be written
{

	int i,j;
	prepareFile(folder, "/QLAREASlope.dat");
	output1 << "L(mm)"<<",    " << "LL(mm)"<<",    " <<"Slope(_)" << ",    " << "Q(mm3/h)" <<",   "<< "AREA(mm2)" <<",   "<< "power(sqrt(Q)*slope)" << endl;
	for(i = 0; i<percentage*M*N; i++)              // it writes a percentage of points: x(%)*MN
	{
		j = int((double)rand()/(double)RAND_MAX*M*N);
		
		output1 << L[j]<<",    "<< LL[j]<<",    " <<SLOPE[j] << ",    " << Q[j] <<",   "<< AREA[j]  <<",   "<<sqrt(Q[j])*SLOPE[j]<< endl;
	}


	writeExceedance(Q, PHI, 101, folder, "/ExceedanceQ.dat"); //it writes the exceedance probabilities of a specific parameter
	writeExceedance(L, PHI, 101, folder, "/ExceedanceL.dat");
	writeExceedance(CURVE, PHI, 101, folder, "/ExceedanceCURVE.dat");
	writeExceedance(SLOPE, PHI, 101, folder, "/ExceedanceSLOPE.dat");
	writeExceedance(AREA, PHI, 101, folder, "/ExceedanceAREA.dat");
	
}

void functionSet::writeExceedance(double *paramMain, double *paramTemp, int NumPoints, string folderName, string fileName) //it writes the exceedance probabilities of a specific parameter
{
	int i,j;
	double max, min;
	int counter, k;
	prepareFile( folderName, fileName);
	
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
	{
		if(paramMain[j*M + i] > max)
			max = paramMain[j*M + i];
		else if(paramMain[j*M + i] < min)
			min = paramMain[j*M + i];
	}
	
	paramTemp[0] = min;
	output1 << min <<"   "<< 1.0 << "\n";
	for(k = 1; k<NumPoints; k++)
	{
		counter = 0;
		paramTemp[k] = min + (double)rand()/(double)RAND_MAX*(max-min);
	
		for(j = 0; j<N; j++)
		for(i = 0; i<M; i++)
		{
			if(paramMain[j*M + i] >= paramTemp[k])
				counter++;
		}
		output1 << paramTemp[k] <<"  "<< counter/double(M*N) << endl;
	}

} 


void functionSet::prepareFile(string folderName, string fileName) //it prepares the output 1 for writing
{
	output1.close();
	folderName += fileName;
	remove(folderName.c_str());
	output1.open(folderName.c_str(),ios::app);
}




void functionSet::filter33(double *param, double *filtered, int filtTime)//it filters the data by a 3 by 3 cells
{
	int i,j;
	int mm1,nn1,ll;
	double temp;
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
		filtered[j*M + i] = param[j*M + i];

	ll=0;
	while(ll<filtTime)
	{
		for(j = 1; j<N-1; j++)
		for(i = 1; i<M-1; i++)
		{
			temp = 0;
			for(mm1 = -1; mm1<2 ; mm1++)
			for(nn1 = -1; nn1<2 ; nn1++)
			{
				temp += filtered[(j + nn1)*M + i + mm1];
				
			}
			filtered[j*M + i] = temp/9.0;
		}
		ll++;
	}
}
double functionSet::filterCELL(double *param, int i, int j, int n)//it filters param[i,j] by an n by n cells
{
	int mm1,nn1;
	double temp = 0;
	for(mm1 = -n/2; mm1<n/2+1 ; mm1++)
	for(nn1 = -n/2; nn1<n/2+1 ; nn1++)
	{
		temp += param[(j + nn1)*M + i + mm1];
		
	}
	return temp/double(n*n);
	
	
}

double functionSet::filternn(double *param, double *filtered, int n)//it filters param by an n by n cells
{
	
	int i,j;
	int mm1,nn1,ll;
	double temp;
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
		filtered[j*M + i] = param[j*M + i];
	
	for(j = 1; j<N-1; j++)
	{
		if(j<n/2)
			for(i = 1; i<M-1; i++)
				filtered[j*M + i] = filterCELL(filtered, i,j, 2*j + 1);

		else if(j>N-n/2 - 1)
			for(i = 1; i<M-1; i++)
				filtered[j*M + i] = filterCELL(filtered, i,j, 2*(-j + N - 1)+1);
		else
		{
			for(i = 1; i<M-1; i++)
			{
				if(i<n/2)
					filtered[j*M + i] = filterCELL(filtered, i,j, 2*i+1);
				else if(i>M-n/2 - 1)
					filtered[j*M + i] = filterCELL(filtered, i,j, 2*(-i + M - 1)+1);
				else
					filtered[j*M + i] = filterCELL(filtered, i,j, n);
			}
		}

	}

	for(j=0; j<N; j++)
	{
		filtered[j*M + 0] =filtered[j*M + 1];filtered[j*M + M - 1] =filtered[j*M + M - 2];
	}

	for(i=0; i<M; i++)
	{
		filtered[(N-1)*M + i] =filtered[(N-2)*M + i];filtered[0*M + i] =filtered[1*M + i];
	}
	
}

void functionSet::filterNetwork(double *Z1, double *Z2, double *Y, double *Qopt) //it filters the value of Z2(experiment) according to the Z1(model) to have the most match between the model and experiment
{
	int i,j;
	int optFilter;
	double minSQR, optMinSQR; 
	double *Zstar1,*PHI1,*Q1, *Q2;
	Zstar1 = new double[M*N];
	PHI1 = new double[M*N];
	Q1 = new double[M*N];
	Q2 = new double[M*N];
	for(i = 1; i < M*N; i++)
	{
		Zstar1[i] = 0;
		PHI1[i] = 0;
		Q1[i] = 0;
		Q2[i] = 0;
	}
	

	replaceError(Zstar1, Z1);
	
	filling(Zstar1, PHI1, Y);
	direction(Zstar1, PHI1, Y);
	drainageQ(Q1);
	

	replaceError(Zstar1, Z2);
	filling(Zstar1, PHI1, Y);
	direction(Zstar1, PHI1, Y);
	drainageQ(Q2);

	minSQR = 0;
	for(i = 0; i < M*N; i++)
	{
		minSQR += (log10(Q1[i]) - log10(Q2[i]))*(log10(Q1[i]) - log10(Q2[i]));
	}
	optFilter = 0;
	optMinSQR = minSQR;
	replaceError(Qopt, Q2);
	
	
	for(i = 1; i < 200; i++)
	{
		filter33(Z2, Zstar1, i);
		filling(Zstar1, PHI1, Y);
		direction(Zstar1, PHI1, Y);
		drainageQ(Q2);
		minSQR = 0;
		for(j = 0; j < M*N; j++)
		{
			minSQR += (log10(Q1[j]+ 1) - log10(Q2[j]+1))*(log10(Q1[j]+1) - log10(Q2[j]+1));
		}
		cout <<i <<"    " <<minSQR <<endl;
		if(minSQR < optMinSQR)
		{
			optFilter = i;
			optMinSQR = minSQR;
			replaceError(Qopt, Q2);
		}
	}
	cout << "optimum filter = "<< optFilter <<endl;
	//delete[] Zstar1; 
	//delete[] PHI1; 
	//delete[] Q1; 
	//delete[] Q2; 
	
}

void functionSet::findMinMax(double *param, double *min, double *max)
{
	int i,j;
	*min = 1e15;
	*max = -1e15;
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
	{
		if(param[j*M + i] > *max)
			*max = param[j*M + i];
		else if(param[j*M + i] < *min)
			*min = param[j*M + i];
	}

}

void functionSet::normalize(double *X) 
{
	int i,j;
	double ave = average(X);
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
		X[j*M + i] = X[j*M + i]/ave;
}




double functionSet::average(double *X) 
{
	int i,j;
	double ave = 0;
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
		ave += X[j*M + i];
	return(ave/double(M*N));
}

double functionSet::RMS(double *X, int size) //https://brenocon.com/blog/2012/03/cosine-similarity-pearson-correlation-and-ols-coefficients/
{
	int i,j;
	double aveX, RMS;
	aveX = average(X);
	RMS = 0;
	
	for(j = 0; j<size; j++)
	{
		RMS += (X[j] - aveX)*(X[j] - aveX);
	}
	
	return sqrt(RMS);
}

double functionSet::meanSquareError(double *X, double *Y, int size) //https://brenocon.com/blog/2012/03/cosine-similarity-pearson-correlation-and-ols-coefficients/
{
	int i,j;
	double RMS;
	RMS = 0;
	
	for(j = 0; j<size; j++)
	{
		RMS += (X[j] - Y[j])*(X[j] - Y[j]);
	}
	
	return sqrt(RMS/size);
}

double functionSet::correlation(double *X, double *Y, int size) //https://brenocon.com/blog/2012/03/cosine-similarity-pearson-correlation-and-ols-coefficients/
{
	int i,j;
	double aveX, aveY, RMSx, RMSy, corr ;
	aveX = average(X);
	aveY = average(Y);
	RMSx = RMS(X, size);
	RMSy = RMS(Y, size);
	corr = 0;
	for(j = 0; j<size; j++)
	{
		corr += (X[j] - aveX)*(Y[j] - aveY);
	}
	
	
	return corr/(RMSx * RMSy);
}


void functionSet::L_down(double *downL)
{
	int i,j,newPlace;
	for(j=1; j<N-1; j++)	
	for(i=1; i<M-1; i++)
	{
		downL[j*M + i] = 0;
		newPlace = drainTo[j*M + i];
		while (newPlace != -1)
		{

			if(newPlace == 	j*M + i + 1 || newPlace == j*M + i - 1)// the initial lenght of each single cell
				downL[j*M + i] += dx;
			else if(newPlace == (j - 1)*M + i || newPlace == (j - 1)*M + i)
				downL[j*M + i] += dy;
			else
				downL[j*M + i] += sqrt(dx*dx + dy*dy);
			cout<<j*M + i<<"   "<<newPlace<<endl;
			newPlace = drainTo[newPlace];
		}
	}
	cout<<drainTo[newPlace]<<endl;
}



