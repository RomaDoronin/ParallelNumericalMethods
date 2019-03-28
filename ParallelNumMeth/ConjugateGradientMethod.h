#pragma once
#include "CRSMatrix.h"

class ConjugateGradientMethod
{
private:
    double* CulcX(double * x, double alf, double * h);
	double* CulcH(double * r1, double bet, double * h);
	double CulcAlf(double * r, double * h, CRSMatrix & A);
	double* CulcR(double * r, double alf, double * h, CRSMatrix & A);
	double CulcBet(double * r, double * r1);

	int n;

public:
    ConjugateGradientMethod();
    ~ConjugateGradientMethod();

    void Solve(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count);
};

double VectorDifSumVal(double * v1, double * v2, int n)
{
	double sum = 0;

	for (int i = 0; i < n; i++)
	{
		sum += abs(v1[i] - v2[i]);
	}

	return sum;
}

double* VectorSum(double * v1, double * v2, int n)
{}

double* VectorMultConst(double val, double * v1, int n)
{}

double VectorScalarMult(double * v1, double * v2, int n)
{}

double* VectorMultMatrix(CRSMatrix & A, double * v1, int n)
{}