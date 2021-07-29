#include "header.h"



int functionSet::temporalBCset(double *X, double *Y, double *Z,string* filename)
{
	int i,j;
	

	scanTime[0] = 0.25;
	scanTime[1] = 0.50;
	scanTime[2] = 1.0;
	scanTime[3] = 2.0;
	scanTime[4] = 4.0;
	scanTime[5] = 8.0;
	scanTime[6] = 16.0;
	double temp; 
	int mm1, nn1, ll = 0;
	int nn;
	remove("gradients.ods");
	output2.open("gradients.ods");
	double aver, avel;
	for (nn = 0; nn < 7; nn++)
	{
		aver = 0;
		avel = 0;
		readLand(X, Y, Z,filename[nn]);///
		replaceError(param, Z);
		//filter(param, param, 1);
		int mm1, nn1, ll = 0;
		double temp;

		
		

		for(i = 0; i < M ; i++)
		{
			down[nn*M + i] = (param[2*M + i]-param[4*M + i])/(2*dy);
			up[nn*M + i] = (param[(N - 2)*M + i] -param[(N - 4)*M + i])/(2*dy);	
			
		}
	
		for(j = 0; j < N ; j++)
		{
			left[nn*N + j] = (param[j*M + 2] - param[j*M + 4])/(2*dx) ;
			right[nn*N + j] = (param[j*M + M - 3] - param[j*M + M - 5])/(2*dx);	
			output2 << nn <<"   "<<j <<"  "<<left[nn*N + j]<< "  "<< right[nn*N + j]<<endl;
			
		}

		while(ll <10)
		{
			for(i =1; i <  M-1 ; i++)
			{
			down[nn*M + i] = (down[nn*M + i] + down[nn*M + i + 1] + down[nn*M + i - 1])/3.0;
			up[nn*M + i] = (up[nn*M + i] + up[nn*M + i + 1] + up[nn*M + i - 1])/3.0;
				
			}
			for(j =1; j <  N -1 ; j++)
			{
				right[nn*N + j] = (right[nn*N + j] + right[nn*N + j + 1] + right[nn*N + j - 1])/3.0;
			left[nn*N + j] = (left[nn*N + j] + left[nn*N + j + 1] + left[nn*N + j - 1])/3.0;
			}
		ll++;
		}
		//cout << nn <<"  "<<aver/N <<"   "<<avel/N <<endl;
	}
	return 1;
}


int functionSet::dtCalc(double *Q, double *SLOPE)
{
	int i,j;
	int k;
	double max = 0, AS;
	double s[4];
	for(i = 1; i<M-1; i++)
	for(j = 1; j<N-1; j++)
	{
		
		AS = pow(Q[j*M + i], m)*pow(SLOPE[j*M + i], n-1.0)*K[j*M + i]*1.4142/dx+(pow(abs(D[j*M + i + 1] - D[j*M + i - 1])/(2*dx) + 1e-5, 0.6666666667)+pow(abs(D[(j + 1)*M + i] - D[(j - 1)*M + i ])/(2*dy) + 1e-5, 0.6666666667))/pow(dx, 0.6666666667);	
	        if(AS > max)
		{
			max = AS;	
		}
	}
	//cout << "max = "<<max <<endl;
	//dt = CFL*dx/(max*K*1.4142);
	dt = CFL/max;
	//if (K == 0 )
		//dt = 0.01;
	//cout << "dt = "<<dt <<endl;
		
	return 1;
}


int functionSet::PHICalc(double *Q, double *SLOPE, double *PHI)//the covection term + noise trem (deposition) + tectonic
{
	int i,j;
	
	/////////////////////////////////////////////////////////////////////attention
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
	
		PHI[j*M + i] = E - delta/dx*K[j*M + i]*(pow(Q[j*M + i], m)*pow(SLOPE[j*M + i], n) - Omega[j*M + i]); //Note, Omega is calculated in the run loop
		//cout<<Omega[j*M + i]<<"  ";
		if(E - PHI[j*M + i]  < 0 || Q[j*M + i] ==0 || SLOPE[j*M + i] ==0 )//
			PHI[j*M + i] = E;
	}
		
	return 0;
}



int functionSet::RungeKutta(double *Z, double *Zstar, double *Q, double *SLOPE, double *dZstar, double *PHI, double *Y)
{
	int i,j;
	slopeCalc(Z,SLOPE);
	
	PHICalc(Q, SLOPE, PHI);
	
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		dZstar[j*M + i] = dt*PHI[j*M + i];
		Zstar[j*M + i] = Z[j*M + i] + dZstar[j*M + i];
		
	}
	
	slopeCalc(Zstar,SLOPE);
	direction(Zstar, PHI, Y);
	drainageQ(Q);
	PHICalc(Q, SLOPE, PHI);
	
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
		{
			Zstar[j*M + i] = Z[j*M + i] + dZstar[j*M + i]/2.0 + dt*PHI[j*M + i]/2.0;
		}
	
	boundarySet(Zstar, Y);
	
	return 0;
}




int functionSet::crankNicholson(double *Zstar, double *Y)
{
	int i,j;
	double f = 1.0;
	
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
		alpha[j*M + i] = 0.5*dt*D[j*M + i]/(dx*dx);
		
	for(j = 0; j<N; j++)
	{	
		i = 1;
		alpha[j*M + i]*=f;
		i = M - 2;
		alpha[j*M + i]*=f;
	}

	for(i = 0; i<M; i++)
	{	
		j = 1;
		alpha[j*M + i]*=f;
		
		j = N-2;
		alpha[j*M + i]*=f;
	}

	for(j = 1; j<N-1; j++)
	{	
		
		i = 1;
		
		a[j*M + 1] = 0;
		b[j*M + 1] = 1 + alpha[j*M + i]/2.0 + alpha[j*M + i + 1]/2.0;   ///////BC:z = 0 (1 + 2 alpha), BC:dz = 0 (1 + alpha)
		c[j*M + 1] = - alpha[j*M + i]/2 - alpha[j*M + i + 1]/2.0;
		d[j*M + 1] = (alpha[j*M + i]/2.0 + alpha[(j - 1)*M + i]/2.0)*Zstar[(j - 1)*M + i] + (1 - alpha[j*M + i] - alpha[(j - 1)*M + i]/2.0- alpha[(j + 1)*M + i]/2.0)*Zstar[j*M + i]  
			+ (alpha[j*M + i]/2.0 + alpha[(j + 1)*M + i]/2.0)*Zstar[(j + 1)*M + i] - (- alpha[j*M + i]/2 - alpha[j*M + i - 1]/2.0)*leftBC[j] ;	

		i = M - 2;
		a[j*M + M - 2] =- alpha[j*M + i]/2 - alpha[j*M + i - 1]/2.0;
		b[j*M + M - 2] = 1 + alpha[j*M + i]/2.0 + alpha[j*M + i - 1]/2.0;///////BC:z = 0 (1 + 2*alpha), BC:dz = 0 (1 + alpha)
		c[j*M + M - 2] = 0;
		d[j*M + M - 2] = (alpha[j*M + i]/2.0 + alpha[(j - 1)*M + i]/2.0)*Zstar[(j - 1)*M + i] + (1 - alpha[j*M + i] - alpha[(j - 1)*M + i]/2.0- alpha[(j + 1)*M + i]/2.0)*Zstar[j*M + i]  
			+ (alpha[j*M + i]/2.0 + alpha[(j + 1)*M + i]/2.0)*Zstar[(j + 1)*M + i] - (- alpha[j*M + i]/2 - alpha[j*M + i + 1]/2.0)*rightBC[j] ;	
	}

	for(j = 1; j<N-1; j++)
	for(i = 2; i<M-2; i++)
	{

			a[j*M + i] = - alpha[j*M + i]/2 - alpha[j*M + i - 1]/2.0;
			b[j*M + i] = 1 + alpha[j*M + i] +  alpha[j*M + i - 1]/2.0 +  alpha[j*M + i + 1]/2.0 ; //1
			c[j*M + i] = - alpha[j*M + i]/2 - alpha[j*M + i + 1]/2.0;
			d[j*M + i] = (alpha[j*M + i]/2.0 + alpha[(j - 1)*M + i]/2.0)*Zstar[(j - 1)*M + i] + (1 - alpha[j*M + i] - alpha[(j - 1)*M + i]/2.0- alpha[(j + 1)*M + i]/2.0)*Zstar[j*M + i]  
			+ (alpha[j*M + i]/2.0 + alpha[(j + 1)*M + i]/2.0)*Zstar[(j + 1)*M + i];	
	/*
		if R=0
		{
			a[j*M + i] = 0;
			b[j*M + i] = 1; 
			c[j*M + i] = 0;
			d[j*M + i] = Zstar[j*M + i];	
		}
		
		//cout <<"i =" << i << "    " << "j =" << j <<"      "<< "b =" <<b[j*M + i] << endl;*/
	}

	Thomas(a,b,c,d,Zstar,0);
	boundarySet(Zstar, Y);

	
	
	for(i = 1; i<M-1; i++)
	{	
		j = 1;
		a[1*M + i] = 0;
		b[1*M + i] = 1 + alpha[j*M + i]/2.0 + alpha[(j + 1)*M + i]/2.0 ;///////BC:z = 0 (1 + 2*alpha), BC:dz = c (1 + alpha) //the value of c is epsilon or alpha*(mainSlope*dy (it is mentioned in d(i,j)
		c[1*M + i] = - alpha[j*M + i]/2.0 - alpha[(j + 1)*M + i]/2.0;
		d[1*M + i] = (alpha[j*M + i]/2.0 + alpha[j*M + i - 1]/2.0)*Zstar[j*M + i-1] + (1 - alpha[j*M + i] - alpha[j*M + i - 1]/2.0 - alpha[j*M + i + 1]/2.0)*Zstar[j*M + i] +  
                        (alpha[j*M + i]/2.0 + alpha[j*M + i + 1]/2.0)*Zstar[j*M + i+1] - (- alpha[j*M + i]/2.0 - alpha[(j - 1)*M + i]/2.0)*downBC[i];	
		
		j = N-2;
		a[(N-2)*M + i] = - alpha[j*M + i]/2.0 - alpha[(j - 1)*M + i]/2.0;
		b[(N-2)*M + i] = 1 + alpha[j*M + i]/2.0 + alpha[(j - 1)*M + i]/2.0 ; ///////BC:z = 0 (1 + 2*alpha), BC:dz = 0 (1 + alpha)
		c[(N-2)*M + i] = 0;
		d[(N-2)*M + i] = (alpha[j*M + i]/2.0 + alpha[j*M + i - 1]/2.0)*Zstar[j*M + i-1] + (1 - alpha[j*M + i] - alpha[j*M + i - 1]/2.0 - alpha[j*M + i + 1]/2.0)*Zstar[j*M + i] +  
                        (alpha[j*M + i]/2.0 + alpha[j*M + i + 1]/2.0)*Zstar[j*M + i+1] - (- alpha[j*M + i]/2.0 - alpha[(j + 1)*M + i]/2.0)*upBC[i];	
	}
/*
	for(i = M/2-int(25/dx); i <= M/2 + int(25/dx); i++)
	{	
		a[1*M + i] = 0;
		b[1*M + i] = 1 + 2*alpha[j*M + i]; ///////BC:z = 0 (1 + 2*alpha), BC:dz = 0 (1 + alpha)
		c[1*M + i] = - alpha[j*M + i] - alpha[(j + 1)*M + i]/4.0 + alpha[(j - 1)*M + i]/4.0;
		d[1*M + i] = alpha*Zstar[1*M + i - 1] + (1 - 2*alpha)*Zstar[1*M + i] + alpha*Zstar[1*M + i+1] ;

	}
*/
	

	i = M/2;j = 1;
	{	
		a[j*M + i] = 0;
		b[j*M + i] = 1 + alpha[j*M + i] + alpha[(j - 1)*M + i]/2.0 + alpha[(j + 1)*M + i]/2.0; ///////BC:z = 0 (1 + 2*alpha), BC:dz = 0 (1 + alpha)
		c[j*M + i] = - alpha[j*M + i]/2.0 - alpha[(j + 1)*M + i]/2.0;;
		d[j*M + i] = (alpha[j*M + i]/2.0 + alpha[j*M + i - 1]/2.0)*Zstar[j*M + i-1] + (1 - alpha[j*M + i] - alpha[j*M + i - 1]/2.0 - alpha[j*M + i + 1]/2.0)*Zstar[j*M + i] +  
                        (alpha[j*M + i]/2.0 + alpha[j*M + i + 1]/2.0)*Zstar[j*M + i+1];	

	}

	for(i = 1; i<M-1; i++)
	for(j = 2; j<N-2; j++)
	{

		
			a[j*M + i] = - alpha[j*M + i]/2.0 - alpha[(j - 1)*M + i]/2.0;
			b[j*M + i] = 1 + alpha[j*M + i] + alpha[(j - 1)*M + i]/2.0 + alpha[(j + 1)*M + i]/2.0;
			c[j*M + i] = - alpha[j*M + i]/2.0 - alpha[(j + 1)*M + i]/2.0;
			d[j*M + i] = (alpha[j*M + i]/2.0 + alpha[j*M + i - 1]/2.0)*Zstar[j*M + i-1] + (1 - alpha[j*M + i] - alpha[j*M + i - 1]/2.0 - alpha[j*M + i + 1]/2.0)*Zstar[j*M + i] +  
                        (alpha[j*M + i]/2.0 + alpha[j*M + i + 1]/2.0)*Zstar[j*M + i+1];	
		
		
	/*	else
		{
			
			a[j*M + i] = 0;
			b[j*M + i] = 1; 
			c[j*M + i] = 0;
			d[j*M + i] = Zstar[j*M + i];	
		}*/
		
	}

	Thomas(a,b,c,d,Zstar,1);
	boundarySet(Zstar, Y);
		
	return 0;
}


int functionSet::Thomas(double *a,double *b,double *c,double *d,double *Z, int rowCol)
{
	int i,j;
	if(rowCol==0)
	{
		for(j = 1; j<N-1; j++)
		{	
			c[j*M + 1] = c[j*M + 1]/b[j*M + 1];
			d[j*M + 1] = d[j*M + 1]/b[j*M + 1];
			b[j*M + 1] = 1.0;

			for(i = 2; i <= M-2; i++)
			{	
				c[j*M + i] = c[j*M + i]/(b[j*M + i] - a[j*M + i]*c[j*M + i - 1]);
				d[j*M + i] = (d[j*M + i] - a[j*M + i]*d[j*M + i - 1])/(b[j*M + i] - a[j*M + i]*c[j*M + i - 1]);
				b[j*M + i] = 1.0;	
			}

			Z[j*M + M - 2] = d[j*M + M - 2];
			for(i = M - 3; i >= 1; i--)
				Z[j*M + i] = d[j*M + i] - c[j*M + i]*Z[j*M + i + 1];
		}
	
	}

	if(rowCol==1)
	{
		for(i = 1; i<M-1; i++)
		{	
			c[1*M + i] = c[1*M + i]/b[1*M + i];
			d[1*M + i] = d[1*M + i]/b[1*M + i];
			b[1*M + i] = 1.0;
			for(j = 2; j <= N-2; j++)
			{	
				c[j*M + i] = c[j*M + i]/(b[j*M + i] - a[j*M + i]*c[(j - 1)*M + i]);
				d[j*M + i] = (d[j*M + i] - a[j*M + i]*d[(j - 1)*M + i])/(b[j*M + i] - a[j*M + i]*c[(j - 1)*M + i]);
				b[j*M + i] = 1.0;
			}

			Z[(N - 2)*M + i] = d[(N - 2)*M + i];
			
			for(j = N - 3; j >=1 ; j--)
				Z[j*M + i] = d[j*M + i] - c[j*M + i]*Z[(j + 1)*M + i];
		}
	}
	return 0;
}

// normal_distribution example



int functionSet::addNoise(double *Z,double dt)
{
	double random,rand1;
	double noise;
	int i,j;
	std::random_device rd;
    	std::default_random_engine generator(rd());
  	std::cauchy_distribution<double> cauchy(0,0.11757900431008157);
	for(i = 1; i<M-1; i++)                                             //adding the noise into the field
	for(j = 1; j<N-1; j++)
	{
		random = 10;
		while (random>1)
			random = abs(cauchy(generator));
		noise = noiseAve*random*dt;
		rand1 = (double) rand() / (RAND_MAX); 
		if(rand1>0.5)
			Z[j*M + i] = Z[j*M + i] + noise;
		else if(rand1<0.5)
			Z[j*M + i] = Z[j*M + i] - noise;
	}
	return 0;
}


int functionSet::boundarySet(double *Z, double *Y)
{
	int i,j;
	
	double timeGap = scanTime[BCstep + 1] - scanTime[BCstep];
	double temp;
	
	for(i = 1; i < M - 1; i++)
	{
		temp = down[BCstep*M + i]+ (down[(BCstep + 1)*M + i] - down[BCstep*M + i])/timeGap*(t - scanTime[BCstep]);
		if(temp > 0 )
			downBC[i] =  temp*dy; //attention downBC = o
		else
			downBC[i] = 0;

		
		temp = Z[1*M + i] + downBC[i];

		if(temp <= Z[0*M + i])
			Z[0*M + i] = temp;
		else if (Z[0*M + i] >Z[1*M + i])
		downBC[i] = Z[0*M + i] - Z[1*M + i];
		else
			downBC[i] = 0;



		temp = up[BCstep*M + i] + (up[(BCstep + 1)*M + i] - up[BCstep*M + i])/timeGap*(t - scanTime[BCstep]);
		if(temp > 0 )
			upBC[i] = temp*dy;
		else
			upBC[i] = 0;

		temp = Z[(N - 2)*M + i] + upBC[i];

		if(temp <= Z[(N - 1)*M + i])
			Z[(N - 1)*M + i] = temp;
		else if (Z[(N - 1)*M + i] >Z[(N - 2)*M + i])
			upBC[i] = Z[(N - 1)*M + i] - Z[(N - 2)*M + i];
		else
			upBC[i] = 0;
	


		/*if(D[j*M + M - 1] + D[j*M + M - 2] >1e-6)
		Z[0*M + i] = Z[1*M + i] + gamma*(D[1*M + i] + D[2*M + i])/(D[0*M + i] + D[1*M + i])*(Z[1*M + i] - Z[2*M + i]); 
		if(Z[0*M + i]< Z[1*M + i])
			Z[0*M + i]= Z[1*M + i]  ;
		
		if(D[(N - 1)*M + i] + D[(N - 2)*M + i] >1e-6)
		Z[(N - 1)*M + i] = Z[(N - 2)*M + i] + gamma*(D[(N - 2)*M + i] + D[(N - 3)*M + i])/(D[(N - 1)*M + i] + D[(N - 2)*M + i])*(Z[(N - 2)*M + i] - Z[(N - 3)*M + i]); 
		if(Z[(N - 1)*M + i]< Z[(N - 2)*M + i] )
			Z[(N - 1)*M + i]= Z[(N - 2)*M + i];*/
		
	}

	for(i = M/2-int(200/dx); i <M/2+int(200/dx); i++)
	{

	downBC[i] = 0.0;
		
		Z[0*M + i] = Z[1*M + i] + downBC[i];
	}
	
		
	
	for(j = 1; j < N - 1 ; j++)
	{

		temp = left[BCstep*N + j] + (left[(BCstep + 1)*N + j] - left[BCstep*N + j])/timeGap*(t - scanTime[BCstep]);
		if(temp > 0 )
			leftBC[j] = temp*dx;
		else
			leftBC[j] = 0;

		temp = Z[j*M + 1] + leftBC[j];

		if(temp <= Z[j*M + 0])
			Z[j*M + 0] = temp;
		else if (Z[j*M + 0] >Z[j*M + 1])
			leftBC[i] = Z[j*M + 0] - Z[j*M + 1];
		else
			leftBC[i] = 0;



		temp = right[BCstep*N + j] + (right[(BCstep + 1)*N + j] - right[BCstep*N + j])/timeGap*(t - scanTime[BCstep]);
		if(temp > 0 )
			rightBC[j] = temp*dx;
		else
			rightBC[j] = 0;

		temp = Z[j*M + M - 2] + rightBC[j];

		if(temp <= Z[j*M + M - 1])
			Z[j*M + M - 1] = temp;
		else if (Z[j*M + M - 1] >Z[j*M + M - 2])
			leftBC[i] = Z[j*M + M - 1] - Z[j*M + M - 2];
		else
			leftBC[i] = 0;

		//rightBC[j] = 0.14;

		/*if(D[j*M + 0] + D[j*M + 1]>1e-6 )
		Z[j*M + 0] = Z[j*M + 1] + gamma*(D[j*M + 1] + D[j*M + 2])/(D[j*M + 0] + D[j*M + 1])*(Z[j*M + 1] - Z[j*M + 2]); 
		if(Z[j*M + 0]< Z[j*M + 1] )
			Z[j*M + 0]= Z[j*M + 1];

		if(D[j*M + M - 1] + D[j*M + M - 2] >1e-6)
		Z[j*M + M - 1] = Z[j*M + M - 2] + gamma*(D[j*M + M - 2] + D[j*M + M - 3])/(D[j*M + M - 1] + D[j*M + M - 2])*(Z[j*M + M - 2] - Z[j*M + M - 3]);  

		if(Z[j*M + M - 1]< Z[j*M + M - 2] )
			Z[j*M + M - 1]= Z[j*M + M - 2] ;*/
	}
	Z[0*M + 0] = (Z[1*M + 0] + Z[0*M + 1])/2;
	Z[0*M + M - 1] = (Z[1*M + M - 1] + Z[0*M + M - 2])/2;
	Z[(N - 1)*M + 0] = (Z[(N - 1)*M + 1] + Z[(N - 2)*M + 0])/2;
	Z[(N - 1)*M + M - 1] = (Z[(N - 2)*M +  M - 1] + Z[(N - 1)*M +  M - 2])/2;

	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	 	if (Z[j*M + i] < mainSlope*Y[j*M + i]-25/cos(atan(mainSlope)))
			Z[j*M + i] = mainSlope*Y[j*M + i]-25/cos(atan(mainSlope));

	Z[0*M + M/2] = 0;
	/*for(i = M/2-int(25/dx); i <= M/2 + int(25/dx); i++)
		Z[0*M + i] = 0;

//FOR analytical verifiction
	/*for(i = 1; i < M; i++)
	{
		Z[0*M + i] = 0;
		Z[(N - 1)*M + i] = 0;
	}
	
	for(j = 1; j < N ; j++)
	{
		Z[j*M + 0] = 0;
		Z[j*M + M - 1] = 0;
	}
*/
	return 0;
}

int functionSet::replaceError(double *Z, double *Zstar)//Z* to Z
{
	int i,j;
	error = 0;
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		//error += (Zstar[j*M + i] - Z[j*M + i])*(Zstar[j*M + i] - Z[j*M + i])/double(M*N);
		Z[j*M + i] = Zstar[j*M + i];
	}
	//cout << error <<endl;
	//error = sqrt(error);
	return 0;
}



	

int functionSet::writeDATA(double *X, double *Y, double *Z, double *X1, double *X2, double *X3, double *X4, double *X5,double *X6,double *X7, string folder, int step)
{
	int i,j;
	char ch[5];
 	double L = 3;
	d2s(double(step),ch,L);
	string str1, str2, str3; 	
	str1 = folder;		
	str1 += "/DATA";
	str1 += ch;
	str1 += ".csv";
	remove(str1.c_str());
	output1.open(str1.c_str(),ios::app);
	double energy = 0;

	/*str2 = folder;	
	str2 += "/SLOPE";
	str2 += ch;
	str2 += ".csv";
	remove(str2.c_str());
	output2.open(str2.c_str(),ios::app);

	str3 = folder;	
	str3 += "/Q";
	str3 += ch;
	str3 += ".csv";
	remove(str3.c_str());
	output3.open(str3.c_str(),ios::app);
	*/
	  

	cout <<"-----------------------WRITING TO:  " <<folder<<"----------------------------------"<<endl;
	//output1<<" X,    Y,    Z,    curve,	filtered_1, filtered_2,	filtered_4, filtered_8, filtered_16" <<endl;
	output1<<" X,    Y,    Z,    H,	slope,Q, AREA,	L, curve,power, downL" <<endl;
			
	//output2<<" X,    Y,    Z,    power" <<endl;
	//output3<<" X,    Y,    Z,    Q" <<endl;
	
	//double random;
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		//random = 0 + (double)rand()/(double)RAND_MAX*0;
		output1 << X[j*M + i] << ",    " << Y[j*M + i]  <<",   " << Z[j*M + i] <<",   " <<X1[j*M + i] <<",   " <<X2[j*M + i]<<",   " <<X3[j*M + i]<<",   " <<X4[j*M + i] <<",   " <<X5[j*M + i] <<",   " <<X6[j*M + i] <<",   " <<X2[j*M + i]*X3[j*M + i]  <<",   "<<X7[j*M + i] << endl;
		
	}
	
/*	for(i = M/2; i<M; i++)
	for(j = 0; j<N; j++)
	{
		random = 0 + (double)rand()/(double)RAND_MAX*0;
		//cout << random <<endl;
		output1 << X[j*M + i] << ",    " << Y[j*M + i] <<",   " << Z[j*M - i] + random <<",   " <<Z[j*M - i] + random <<",   " <<SLOPE[j*M - i] <<",   " <<Q[j*M - i] <<",   " <<CURVE[j*M - i] <<",   " <<sqrt(Q[j*M - i])*SLOPE[j*M - i] << "   " << endl;
		
	}*/

	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		if(DIR[j*M + i] == 1 || DIR[j*M + i] == 4 || DIR[j*M + i] == 16 || DIR[j*M + i] == 32)
		{
	 		energy += sqrt(X3[j*M + i])*dx;
		}
		else
		{
	 		energy += sqrt(X3[j*M + i]*2.0)*dx;
		}
				
	}
	cout << "energy (Q) = "<<energy<<endl;


	energy = 0;
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		if(DIR[j*M + i] == 1 || DIR[j*M + i] == 4 || DIR[j*M + i] == 16 || DIR[j*M + i] == 32)
		{
	 		energy += sqrt(X4[j*M + i])*dx;
		}
		else
		{
	 		energy += sqrt(X4[j*M + i]*2.0)*dx;
		}
				
	}
	
	cout << "energy (A) = "<<energy<<endl;
	output1.close();
	//output2.close();
	//output3.close();
	return 0;
}


int functionSet::writeXYZ(double *X, double *Y, double *Z, string folder, int step)
{

	int i,j;
	char ch[5];
 	double L = 3;
	d2s(double(step),ch,L);
	string str1, str2, str3; 	
	str1 = folder;		
	str1 += "/XYZ";
	str1 += ch;
	str1 += ".dat";
	remove(str1.c_str());
	output1.open(str1.c_str(),ios::app);
	for(j = 0; j<N; j++)
	for(i = 0; i<M; i++)
	{
		output1 << X[j*M + i] << "   " << Y[j*M + i]  <<"   " << Z[j*M + i] <<endl;	
	}

	output1.close();

}



int functionSet::analytic(double *X, double *Y, double *Z)
{
	double numb = 200;
	int k, i,j;
	remove("analytics.csv");
	output1.open("analytics.csv",ios::app);
	output1<<" X,    Y,    Z,    H" <<endl;
	for(i = 0; i<M; i++)
	for(j = 0; j<N; j++)
	{
		Z[j*M + i] = (1-X[j*M + i]*X[j*M + i])/2.0;
		for(k = 1; k<numb; k = k+2)
			Z[j*M + i] += -16/(PI*PI*PI)*( sin(k*PI*(1+X[j*M + i])/2)/(k*k*k*sinh(k*PI))*(sinh(k*PI*(1+Y[j*M + i])/2) + sinh(k*PI*(1-Y[j*M + i])/2)));
			
		output1 << X[j*M + i] << ",    " << Y[j*M + i] <<",   " << Z[j*M + i] <<",   " <<Z[j*M + i] << "   " << endl;	
	}
	output1.close();
	return 0;
}




void functionSet::d2s(double num, char s[], unsigned int width)

{
	int i,j;
	//--width; // allow room for \0.
	double decimal = abs(num) - abs(int(num));
	int integer = abs(int(num));
	deque<char> result;
         if (integer == 0) 
		{result.push_front('0');}
	while (integer != 0)
	{
		result.push_front(integer % 10 + '0');
		integer /= 10;
	}
       

	if (num < 0) // check for negative value
		result.push_front('-');
	//if (decimal != 0) // check for decimal value
		result.push_back('.');
		//if(int(decimal*10) == 0.0){result.push_back(0 + '0'); decimal *= 10;width +=2;}

	while (decimal != 0.00 && result.size() < width)
	{
		decimal *= 10.0;
		
		result.push_back(int(decimal) % 10 + '0');
		decimal = decimal - int(decimal);
	}
	while (result.size() < width)
		result.push_front('0'); // pad with spaces.
	for (int i=0; i<width; ++i) // results into array s
		s[i] = result[i];
	s[width] = '\0';
}




