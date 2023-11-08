/*
 *  MatrixMathExt.h Library for Matrix Math Extended
 *
 *  Created by AranaCorp on 18/06/2017.
 *  Modified from Charlie Matlack's MatrixMath library
 *  summary:
 *		- modified PseudoInvert output
 *		- add PseudoDet
 *		- add Determinant
 *		- add crossProduct
 */

#ifndef MatrixMathExt_h
#define MatrixMathExt_h

#if defined(ARDUINO) && ARDUINO >= 100
#include "Arduino.h"
#else
#include "WProgram.h"
#endif

class MatrixMathExt
{
public:
	//MatrixMathExt();
	void Print(float* A, int m, int n, String label);
	void Copy(float* A, int n, int m, float* B);
	void Multiply(float* A, float* B, int m, int p, int n, float* C);
	void Add(float* A, float* B, int m, int n, float* C);
	void Subtract(float* A, float* B, int m, int n, float* C);
	void Transpose(float* A, int m, int n, float* C);
	void Scale(float* A, int m, int n, float k);
	int Invert(float* A, int n);
	void PseudoInvert(float* J, int m,int n, float* invJ);
	double PseudoDet(float* J,int m, int n);
	double Determinant(float* A, int n);
	void CrossProduct1D(float* a, float* b, float* r);
private:
	double CalcDeterminant( float **mat, int order);	
	int GetMinor(float **src, float **dest, int row, int col, int order);
};

extern MatrixMathExt Matrix;
#endif
