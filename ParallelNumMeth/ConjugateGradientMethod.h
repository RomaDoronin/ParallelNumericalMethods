#pragma once
#include "CRSMatrix.h"
#include <omp.h>

class ConjugateGradientMethod
{
private:
    void CulcX(double * x, double alf, double * h, double * culcRes);
    void CulcH(double * r1, double bet, double * h, double * culcRes);
    double CulcAlf(double * r, double * h, CRSMatrix & A);
    void CulcR(double * r, double alf, double * h, CRSMatrix & A, double * culcRes);
    double CulcBet(double * r, double * r1);

    int n;
	int ProcNum;
	
	double * vectorMultMatrix;
	double * vectorMultConst;

public:
    ConjugateGradientMethod(int _ProcNum);
    ~ConjugateGradientMethod();

    void Solve(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count);
};

double VectorDifSumVal(double * v1, double * v2, int n, int ProcNum);

void VectorSum(double * v1, double * v2, int n, double * vRes);

void VectorMultConst(double val, double * v1, int n, double * vRes);

double VectorScalarMult(double * v1, double * v2, int n, int ProcNum);

void VectorMultMatrix(CRSMatrix & A, double * v1, int n, double * vRes);

void Vec1CopyToVec2(double * v1, double * v2, int n);
