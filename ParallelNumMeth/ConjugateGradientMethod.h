#pragma once
#include "CRSMatrix.h"
#include <omp.h>

class ConjugateGradientMethod
{
private:
    void CulcX(double * x, double alf, double * h, double * culcRes);
    void CulcH(double * r1, double bet, double * h, double * culcRes);
    double CulcAlf(double * r, double * h, const CRSMatrix & A);
    void CulcR(double * r, double alf, double * h, const CRSMatrix & A, double * culcRes);
    double CulcBet(double * r, double * r1);

    int n;
	int ProcNum;
	
	double * vectorMultMatrix;
	double * vectorMultConst;

public:
    ConjugateGradientMethod(int _ProcNum);
    ~ConjugateGradientMethod();

    void Solve(const CRSMatrix & A, const double * b, double eps, int max_iter, double * x, int & count);
};

double VectorDifSumVal(const double * v1, const double * v2, int n, int ProcNum);
double VectorDifSumVal(const double * v1, const double * v2, int n);

void VectorSum(const double * v1, const double * v2, int n, double * vRes);

void VectorMultConst(double val, const double * v1, int n, double * vRes);

double VectorScalarMult(double * v1, double * v2, int n, int ProcNum);

void VectorMultMatrix(const CRSMatrix & A, const double * v1, int n, double * vRes);

void Vec1CopyToVec2(const double * v1, double * v2, int n);
