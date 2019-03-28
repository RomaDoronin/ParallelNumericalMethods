#include "ConjugateGradientMethod.h"


double* ConjugateGradientMethod::CulcX(double * x, double alf, double * h)
{
	return VectorSum(x, VectorMultConst(alf, h, n), n);
}

double* ConjugateGradientMethod::CulcH(double * r1, double bet, double * h)
{
	return VectorSum(r1, VectorMultConst(bet, h, n), n);
}

double ConjugateGradientMethod::CulcAlf(double * r, double * h, CRSMatrix & A)
{
	return VectorScalarMult(r, r, n) / VectorScalarMult(VectorMultMatrix(A, h, n), h, n);
}

double* ConjugateGradientMethod::CulcR(double * r, double alf, double * h, CRSMatrix & A)
{
	return VectorSum(r, VectorMultConst(-alf, VectorMultMatrix(A, h, n), n), n);
}

double ConjugateGradientMethod::CulcBet(double * r, double * r1)
{
	return VectorScalarMult(r1, r1, n) / VectorScalarMult(r, r, n);
}

ConjugateGradientMethod::ConjugateGradientMethod()
{
}

ConjugateGradientMethod::~ConjugateGradientMethod()
{
}

void ConjugateGradientMethod::Solve(CRSMatrix & A, double * b, double eps, int max_iter, double * x, int & count)
{
	n = A.GetN();

	double * r;
	double * rPr;
	double * h;
	double alf;
	double bet;
	double * xPr;

	// Начальные данные
	r = VectorSum(b, VectorMultConst(-1, VectorMultMatrix(A, x, n), n), n);
	h = r;
	alf = VectorScalarMult(r, r, n) / VectorScalarMult(VectorMultMatrix(A, h, n), h, n);
	xPr = x;
	x = CulcX(x, alf, h);
	count = 1;
	
	while (VectorDifSumVal(x, xPr, n) > eps || max_iter >= count)
	{
		rPr = r;
		r = CulcR(rPr, alf, h, A); // R_s
		bet = CulcBet(rPr, r);     // BET_s-1
		h = CulcH(r, bet, h);      // H_s
		alf = CulcAlf(r, h, A);    // ALF_s
		xPr = x;
		x = CulcX(x, alf, h);      // X_s+1

		count++;
	}
}
