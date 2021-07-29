#include "header.h"


void functionSet::Ktest(double *X, double *Y, double *center, int maxK, string folderName, string fileName) 
{
	
	int i, j, k, myK, n, m;
	int *group, *groupPop;
	group = new int [M*N]; 
	groupPop = new int [maxK]; 

	double error, min, distance, aa;
	double *bb;
	bb = new double [maxK]; 
	
	
	output2.close();
	folderName += fileName;
	remove(folderName.c_str());
	output2.open(folderName.c_str(),ios::app);
	
	for(myK=2; myK<maxK; myK+=1)
	{
		error = 0;
			
		Kmean(X, Y, center,M*N, myK, folderName, "/test.txt");
		
		for(k = 0; k < myK; k++)
		{
			groupPop[k] = 0;

		}
		
		for(j = 0; j<N; j++)
		for(i = 0; i<M; i++)
		{
			min = 1e15;
		  	for(k = 0; k < myK; k++)
		  	{
			   	//distance = sqrt( (X[j*M + i] - center[0*myK + k])*(X[j*M + i] - center[0*myK + k]) + (Y[j*M + i] - center[1*myK + k])*(Y[j*M + i] - center[1*myK + k]) );
				distance = sqrt( (X[j*M + i] - center[0*myK + k])*(X[j*M + i] - center[0*myK + k]) )/1e15;
				
				if(distance < min)
				{
					group[j*M + i] = k;
					min = distance;
				}
	 	 	}
			for(k = 0; k < myK; k++)
			{
				if(group[j*M + i] == k)
					groupPop[k]+=1;

			}
			//cout << group[j*M + i] << endl;
		}
		
		for(j = 0; j<N; j++)
		for(i = 0; i<M; i++)
		{
			for(k = 0; k < myK; k++)
			if(group[j*M + i] == k)
			error += sqrt( (X[j*M + i] - center[0*myK + k])*(X[j*M + i] - center[0*myK + k]) );
			//error += sqrt( (X[j*M + i] - center[0*myK + k])*(X[j*M + i] - center[0*myK + k]) + (Y[j*M + i] - center[1*myK + k])*(Y[j*M + i] - center[1*myK + k]) );	
		}
		
		output2 << myK <<"   "<< error/(M*N) <<"  ";

		
		error = 0; //to save the summation of S's
		for(j = 0; j<N; j++)
		for(i = 0; i<M; i++)
		{
			a = 0;
			for(k = 0; k < myK; k++)
				bb[k] = 0;

			for(n = 0; n<N; n++) // the distance from each point is added to b[k]'s
			for(m = 0; m<M; m++)
			{
				//bb[group[n*M + m]] += sqrt( (X[j*M + i] - X[n*M + m])*(X[j*M + i] - X[n*M + m]) + (Y[j*M + i] - Y[n*M + m])*(Y[j*M + i] - Y[n*M + m]) );
				bb[group[n*M + m]] += sqrt( (X[j*M + i] - X[n*M + m])*(X[j*M + i] - X[n*M + m]) );

			}
			for(k = 0; k < myK; k++)
				bb[k] /= groupPop[k];
			
			aa = bb[group[j*M + i]];
			
			min = 1e20;
			for(k = 0; k < myK && k != group[j*M + i]; k++)
				if(bb[k] < min)
				{
					min = bb[k];
				}


			if(aa < min)
				error += 1- aa/min;
			else if ( aa > min)
				error += min/aa - 1;
			
				
		}
		
		output2 << myK <<"   "<< error/(M*N) <<endl;
 //the average S

	}
	
	
	delete[] group;
	delete[] bb;
	delete[] groupPop;
}


void functionSet::Kmean(double *XX, double *YY, double *center, int NumPoint, int myK, string folderName, string fileName) //Attention:the distance can be changed in X, Y or both directions, also be careful abou the calculation of variance (just variance of y value has been calculated.
{
	int i,j;
	double distance, *groupVar, *centerOLD, *X,*Y, min,max;
	int *group, *groupPop,k;
	
	group = new int [NumPoint]; 
	X = new double [NumPoint]; 
	Y = new double [NumPoint]; 
	centerOLD = new double [2*myK]; 
	groupPop = new int [myK];
	groupVar = new double [myK];


	for(j = 0; j<NumPoint; j++)
	{
		X[j] = log(XX[j])/log(10.0);//log(XX[j])/log(10.0);
		Y[j] = log(YY[j])/log(10.0);//log(YY[j])/log(10.0);
	}
	
	min = 1e15;
	max = -1e15;
	for(j = 0; j<NumPoint; j++)
	{
		if(X[j] > max)
			max = X[j];
		else if(X[j] < min)
			min = X[j];
	}

	

	/*	for(j = 0; j<N; j++)
		for(i = 0; i<M; i++)
		{
			if(Y[j*M + i] <= 1e-15)
				Y[j*M + i] = 1e-15;
			if(X[j*M + i] <= 1e-15)
				X[j*M + i] = 1e-15;
			X[j*M + i] = log(X[j*M + i])/log(10.0);
			Y[j*M + i] = log(Y[j*M + i])/log(10.0);
			cout << X[j*M + i] << "   "<<Y[j*M + i] << endl;
		}*/
				

 	
		

	 	for(k=0; k<myK; k++)// initial group centers
		{
			center[0*myK + k] = min + (max-min)/myK*k;	//uniform distribution along X
			i = int((double)rand()/(double)RAND_MAX*NumPoint);
			center[1*myK + k] = Y[i];	
		}
		
	
		
	 double ok = 0; 
	 while(ok < myK - 1)
	 {
		ok = 0;
		for(j = 0; j<NumPoint; j++)
		{
			min = 1e15;
		  	for(k = 0; k < myK; k++)
		  	{
			   	//distance = abs( (X[j*M + i] - center[0*myK + k])*(X[j*M + i] - center[0*myK + k]) + (Y[j*M + i] - center[1*myK + k])*(Y[j*M + i] - center[1*myK + k]) )/1e15;
				distance = sqrt( (X[j] - center[0*myK + k])*(X[j] - center[0*myK + k]) )/1e15;
				
				if(distance < min)
				{
					group[j] = k;
					min = distance;
				}
	 	 	}
			//cout << group[j*M + i] << endl;
		}

		for(k = 0; k < myK; k++)
		{
			centerOLD[0*myK + k] = center[0*myK + k];
			centerOLD[1*myK + k] = center[1*myK + k];
			center[0*myK + k] = 0;
			center[1*myK + k] = 0;
			groupPop[k] = 0;
			for(j = 0; j<NumPoint; j++)
			{
				if(group[j] == k)
				{
					center[0*myK + k] += X[j];	
					center[1*myK + k] += Y[j];
					groupPop[k]+=1;	
				}
			}
			
			if(groupPop[k] >= 1)
			{
				center[0*myK + k]/= double(groupPop[k]);
				center[1*myK + k]/= double(groupPop[k]);	
			}
			else
			{
				i = int((double)rand()/(double)RAND_MAX*M*N);
				center[0*myK + k] = X[i];	
				center[1*myK + k] = Y[i];
			}
			
			
			//cout << groupPop<<"   "<<center[0*myK + k] <<"   " << center[1*myK + k] << endl;
		}
		
		for(k = 0; k < myK; k++)
			if (centerOLD[1*k] == center[0*myK + k] || centerOLD[2*k] == center[1*myK + k])
			ok++;
		cout << ok <<endl;	
				
	}

	for(k = 0; k < myK; k++)
	{
		groupVar[k] = 0;
		for(j = 0; j<NumPoint; j++)
		{
			if(group[j] == k)
				groupVar[k] += (center[1*myK + k] - Y[j])*(center[1*myK + k] - Y[j])/double(groupPop[k]);	
			
		}

		

	}
		  	

	prepareFile(folderName, fileName);
	output1 << "Xmin" <<",   " << "Ymin"  <<",      "<< "groupPop"<<",   "<< "Variance"<< endl;
	for(int k = 0; k < myK; k++)
		output1 << pow(10, center[0*myK + k]) <<",   " <<pow(10.0,center[1*myK + k]) <<",   " << groupPop[k]<< ",   " << groupVar[k] << endl;
	
	delete[] group;
	delete[] groupPop;
	delete[] centerOLD;
	delete[] groupVar;

}









