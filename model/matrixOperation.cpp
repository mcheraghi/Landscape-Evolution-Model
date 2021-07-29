#include "header.h"



int matrixOperations::transpose(double *A, double *B, int row, int col)
{
		
   	int i,j;
	double *temp;
	temp = new double[row*col];
	

	for(i=0; i<row; i++)
	for(j=0; j<col; j++)
	{
		temp[i*col + j]  = A[j*row + i];
	}
	
	for(i=0; i<row; i++)
	for(j=0; j<col; j++)
	{
		B[i*col + j]  	= temp[i*col + j];
	}

	delete [] temp;

	
	return 0;
}

int matrixOperations::multiple(double *A, double *B, double *C, int Ar, int Ac, int Br, int Bc)
{
		
   	int i,j,k;
/* If colum of first matrix in not equal to row of second matrix, asking user to enter the size of matrix again. */
	while (Ac!=Br)
    	{
		cout << "Error! column of first matrix not equal to row of second.";
		return 1;
    	}

/* Initializing elements of matrix mult to 0.*/
	for(i=0; i<Ar; ++i)
	for(j=0; j<Bc; ++j)
	{
		C[j*Ar + i] = 0;
	}

/* Multiplying matrix a and b and storing in array mult. */
	for(i=0; i<Ar; ++i)
	for(j=0; j<Bc; ++j)
	for(k=0; k<Ac; ++k)
	{
		
		C[j*Ar + i] += A[k*Ar + i]*B[j*Br + k];
	}

	return 0;
}


/// calculate the cofactor of element (row,col)
int matrixOperations::getMinor(double **src, double **dest, int row, int col, int order)
{
    // indicate which col and row is being copied to dest
    int colCount=0,rowCount=0;

    for(int i = 0; i < order; i++ )
    {
        if( i != row )
        {
            colCount = 0;
            for(int j = 0; j < order; j++ )
            {
                // when j is not the element
                if( j != col )
                {
                    dest[rowCount][colCount] = src[i][j];
                    colCount++;
                }
            }
            rowCount++;
        }
    }

    return 1;
}

// Calculate the determinant recursively.
double matrixOperations::calcDeterminant( double **mat, int order)
{
    // order must be >= 0
	// stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];

    // the determinant value
    double det = 0;

    // allocate the cofactor matrix
    double **minor;
    minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = new double[order-1];

    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        getMinor( mat, minor, 0, i , order);
        // the recusion is here!
        det += pow( -1.0, i ) * mat[0][i] * calcDeterminant( minor,order-1 );
    }

    // release memory
    for(int i=0;i<order-1;i++)
        delete [] minor[i];
    delete [] minor;

    return det;
}

// matrix inversioon
// the result is put in Y
void matrixOperations::matrixInversion(double *A, int order, double *Y)
{
	// memory allocation
    double *temp1 = new double[order*order];
    double **AA = new double*[order];
    for(int i=0;i<order;i++)
        AA[i] = temp1+(i*(order));

	for(int j=0;j<order;j++)
        for(int i=0;i<order;i++)
		AA[i][j] = A[j*order + i];

    // get the determinant of a
	//cout<<calcDeterminant(AA,order)<<endl;
    double det = 1.0/calcDeterminant(AA,order);

    // memory allocation
    double *temp = new double[(order-1)*(order-1)];
    double **minor = new double*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = temp+(i*(order-1));

    for(int j=0;j<order;j++)
    {
        for(int i=0;i<order;i++)
        {
            // get the co-factor (matrix) of A(j,i)
            getMinor(AA,minor,i,j,order);
            Y[j*order + i] = det*calcDeterminant(minor,order-1);
            if( (i+j)%2 == 1)
                Y[j*order + i] = -Y[j*order + i];
        }
    }

    // release memory
    delete [] minor[0];
    delete [] minor;
}

