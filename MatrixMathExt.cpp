/*
 *  MatrixMathExt.cpp Library for Matrix Math Extended
 *
 *  Created by AranaCorp on 18/06/2017.
 *  Modified from Charlie Matlack's MatrixMath library
 *
 */

#include "MatrixMathExt.h"

#define NR_END 1

MatrixMathExt Matrix;			// Pre-instantiate

// Matrix Printing Routine
// Uses tabs to separate numbers under assumption printed float width won't cause problems
void MatrixMathExt::Print(float* A, int m, int n, String label){
	// A = input matrix (m x n)
	int i,j;
	Serial.println();
	Serial.println(label);
	for (i=0; i<m; i++){
		for (j=0;j<n;j++){
			Serial.print(A[n*i+j]);
			Serial.print("\t");
		}
		Serial.println();
	}
}

void MatrixMathExt::Copy(float* A, int n, int m, float* B)
{
	int i, j, k;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			B[n*i+j] = A[n*i+j];
		}
}

//Matrix Multiplication Routine
// C = A*B
void MatrixMathExt::Multiply(float* A, float* B, int m, int p, int n, float* C)
{
	// A = input matrix (m x p)
	// B = input matrix (p x n)
	// m = number of rows in A
	// p = number of columns in A = number of rows in B
	// n = number of columns in B
	// C = output matrix = A*B (m x n)
	int i, j, k;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
		{
			C[n*i+j]=0;
			for (k=0;k<p;k++)
				C[n*i+j]= C[n*i+j]+A[p*i+k]*B[n*k+j];
		}
}


//Matrix Addition Routine
void MatrixMathExt::Add(float* A, float* B, int m, int n, float* C)
{
	// A = input matrix (m x n)
	// B = input matrix (m x n)
	// m = number of rows in A = number of rows in B
	// n = number of columns in A = number of columns in B
	// C = output matrix = A+B (m x n)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			C[n*i+j]=A[n*i+j]+B[n*i+j];
}


//Matrix Subtraction Routine
void MatrixMathExt::Subtract(float* A, float* B, int m, int n, float* C)
{
	// A = input matrix (m x n)
	// B = input matrix (m x n)
	// m = number of rows in A = number of rows in B
	// n = number of columns in A = number of columns in B
	// C = output matrix = A-B (m x n)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			C[n*i+j]=A[n*i+j]-B[n*i+j];
}


//Matrix Transpose Routine
void MatrixMathExt::Transpose(float* A, int m, int n, float* C)
{
	// A = input matrix (m x n)
	// m = number of rows in A
	// n = number of columns in A
	// C = output matrix = the transpose of A (n x m)
	int i, j;
	for (i=0;i<m;i++)
		for(j=0;j<n;j++)
			C[m*j+i]=A[n*i+j];
}

void MatrixMathExt::Scale(float* A, int m, int n, float k)
{
	for (int i=0; i<m; i++)
		for (int j=0; j<n; j++)
			A[n*i+j] = A[n*i+j]*k;
}


//Matrix Inversion Routine
// * This function inverts a matrix based on the Gauss Jordan method.
// * Specifically, it uses partial pivoting to improve numeric stability.
// * The algorithm is drawn from those presented in 
//	 NUMERICAL RECIPES: The Art of Scientific Computing.
// * The function returns 1 on success, 0 on failure.
// * NOTE: The argument is ALSO the result matrix, meaning the input matrix is REPLACED
int MatrixMathExt::Invert(float* A, int n)
{
	// A = input matrix AND result matrix
	// n = number of rows = number of columns in A (n x n)
	int pivrow;		// keeps track of current pivot row
	int k,i,j;		// k: overall index along diagonal; i: row index; j: col index
	int pivrows[n]; // keeps track of rows swaps to undo at end
	float tmp;		// used for finding max value and making column swaps

	for (k = 0; k < n; k++)
	{
		// find pivot row, the row with biggest entry in current column
		tmp = 0;
		for (i = k; i < n; i++)
		{
			if (abs(A[i*n+k]) >= tmp)	// 'Avoid using other functions inside abs()?'
			{
				tmp = abs(A[i*n+k]);
				pivrow = i;
			}
		}

		// check for singular matrix
		if (A[pivrow*n+k] == 0.0f)
		{
			Serial.println("Inversion failed due to singular matrix");
			return 0;
		}

		// Execute pivot (row swap) if needed
		if (pivrow != k)
		{
			// swap row k with pivrow
			for (j = 0; j < n; j++)
			{
				tmp = A[k*n+j];
				A[k*n+j] = A[pivrow*n+j];
				A[pivrow*n+j] = tmp;
			}
		}
		pivrows[k] = pivrow;	// record row swap (even if no swap happened)

		tmp = 1.0f/A[k*n+k];	// invert pivot element
		A[k*n+k] = 1.0f;		// This element of input matrix becomes result matrix

		// Perform row reduction (divide every element by pivot)
		for (j = 0; j < n; j++)
		{
			A[k*n+j] = A[k*n+j]*tmp;
		}

		// Now eliminate all other entries in this column
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				tmp = A[i*n+k];
				A[i*n+k] = 0.0f;  // The other place where in matrix becomes result mat
				for (j = 0; j < n; j++)
				{
					A[i*n+j] = A[i*n+j] - A[k*n+j]*tmp;
				}
			}
		}
	}

	// Done, now need to undo pivot row swaps by doing column swaps in reverse order
	for (k = n-1; k >= 0; k--)
	{
		if (pivrows[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				tmp = A[i*n+k];
				A[i*n+k] = A[i*n+pivrows[k]];
				A[i*n+pivrows[k]] = tmp;
			}
		}
	}
	return 1;
}
//MatrixMath++

//Matrix Pseudo Inversion Routine
// * This function generate a pseudo inverse matrix
void MatrixMathExt::PseudoInvert(float* J, int m,int n, float* invJ){

	float* Jt;
	float* C;
	
    if (m>n){
        //invJ=inv(transpose(J)*J)*transpose(J);
        Transpose(J,m,n,Jt);
        Multiply(Jt, J, n, m, n, C);
        Invert(C,n);
        Multiply(C,Jt,n,n,m,invJ);
    }else if (m<n){
        //invJ=transpose(J)*inv(J*transpose(J));
        Transpose(J,m,n,Jt);
        Multiply(J, Jt, m, n, m, C);
        Invert(C,m);
        Multiply(Jt,C,n,m,m,J);
    }else{
        //invJ=Invert(J);
        Copy(J, m, m, invJ);
        Invert(invJ,m);
    }
}

double MatrixMathExt::PseudoDet(float* J,int m, int n){
    //function detJ=pseudoDet(J)
    //Pseudo determinant
    float* Jt;
	float* C;
    double detJ=0;

    if (m>n){
        //detJ=det(transpose(J)*J);
        Transpose(J,m,n,Jt);
        Multiply(Jt,J, n, m, n, C);
        detJ=Determinant(C, n);
    }else if (m<n){
        //detJ=det(J*transpose(J));
        Transpose(J,m,n,Jt);
        Multiply(J,Jt, m, n, m, C);
        detJ=Determinant(C, m);
    }else{
        //detJ=det(J);
        detJ=Determinant(J, m);
    }
    return detJ;
}

//MathMatrix++

// calculate the cofactor of element (row,col)
int MatrixMathExt::GetMinor(float **src, float **dest, int row, int col, int order){
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
double MatrixMathExt::CalcDeterminant( float **mat, int order){
    // order must be >= 0
    // stop the recursion when matrix is a single element
    if( order == 1 )
        return mat[0][0];

    // the determinant value
    float det = 0;

    // allocate the cofactor matrix
    float **minor;
    minor = new float*[order-1];
    for(int i=0;i<order-1;i++)
        minor[i] = new float[order-1];

    for(int i = 0; i < order; i++ )
    {
        // get minor of element (0,i)
        GetMinor( mat, minor, 0, i , order);
        // the recusion is here!

        det += (i%2==1?-1.0:1.0) * mat[0][i] * CalcDeterminant(minor,order-1);
        //det += pow( -1.0, i ) * mat[0][i] * CalcDeterminant( minor,order-1 );
    }

    // release memory
    for(int i=0;i<order-1;i++)
        delete [] minor[i];
    delete [] minor;

    return det;
}

double MatrixMathExt::Determinant(float* A, int n){
    float **B = new float *[n];
    for(int p=0;p<4;p++)B[p]=new float [n];

    for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) B[i][j]=A[i*n+j];

    return CalcDeterminant(B,n);
    delete [] B;
}

void MatrixMathExt::CrossProduct1D(float* a, float* b, float* r){
    /*Ex: double v1[3]={1,0,0};
    double v2[3]={0,0,1};
    double w[3];
    CrossProduct1D(v1,v2,w);*/

  r[0] = a[1]*b[2]-a[2]*b[1];
  r[1] = a[2]*b[0]-a[0]*b[2];
  r[2] = a[0]*b[1]-a[1]*b[0];

}
