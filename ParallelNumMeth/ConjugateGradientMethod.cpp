#include "ConjugateGradientMethod.h"


void ConjugateGradientMethod::CulcX(double * x, double alf, double * h, double * culcRes)
{
	VectorMultConst(alf, h, n, vectorMultConst);
    VectorSum(x, vectorMultConst, n, culcRes);
}

void ConjugateGradientMethod::CulcH(double * r1, double bet, double * h, double * culcRes)
{
	VectorMultConst(bet, h, n, vectorMultConst);
    return VectorSum(r1, vectorMultConst, n, culcRes);
}

double ConjugateGradientMethod::CulcAlf(double * r, double * h, CRSMatrix & A)
{
	VectorMultMatrix(A, h, n, vectorMultMatrix);
    return VectorScalarMult(r, r, n, ProcNum) / VectorScalarMult(vectorMultMatrix, h, n, ProcNum);
}

void ConjugateGradientMethod::CulcR(double * r, double alf, double * h, CRSMatrix & A, double * culcRes)
{
	VectorMultMatrix(A, h, n, vectorMultMatrix);
	VectorMultConst(-alf, vectorMultMatrix, n, vectorMultConst);
    VectorSum(r, vectorMultConst, n, culcRes);
}

double ConjugateGradientMethod::CulcBet(double * r, double * r1)
{
    return VectorScalarMult(r1, r1, n, ProcNum) / VectorScalarMult(r, r, n, ProcNum);
}

ConjugateGradientMethod::ConjugateGradientMethod(int _ProcNum) : ProcNum(_ProcNum)
{
}

ConjugateGradientMethod::~ConjugateGradientMethod()
{
	delete vectorMultMatrix;
	delete vectorMultConst;
}

void ConjugateGradientMethod::Solve(CRSMatrix & A, double *b, double eps, int max_iter, double *x, int & count)
{
    n = A.GetN();

    double * r = new double[n];
    double * rPr = new double[n];
    double * h = new double[n];
	double * hPr = new double[n];
	double alf;
	double bet;
    double * xPr = new double[n];

    // Начальные данные
	vectorMultMatrix = new double[n];
	VectorMultMatrix(A, x, n, vectorMultMatrix);
	vectorMultConst = new double[n];
	VectorMultConst(-1, vectorMultMatrix, n, vectorMultConst);
    VectorSum(b, vectorMultConst, n, r);
	VectorSum(b, vectorMultConst, n, h);

	VectorMultMatrix(A, h, n, vectorMultMatrix);
    alf = VectorScalarMult(r, r, n, ProcNum) / VectorScalarMult(vectorMultMatrix, h, n, ProcNum);

	Vec1CopyToVec2(x, xPr, n);
    CulcX(xPr, alf, h, x);
    count = 1;
    
    while (VectorDifSumVal(x, xPr, n, ProcNum) > eps && max_iter >= count)
    {
		Vec1CopyToVec2(r, rPr, n);
        CulcR(rPr, alf, h, A, r);  // R_s
        bet = CulcBet(rPr, r);     // BET_s-1
		Vec1CopyToVec2(h, hPr, n);
        CulcH(r, bet, hPr, h);     // H_s
        alf = CulcAlf(r, h, A);    // ALF_s
		Vec1CopyToVec2(x, xPr, n);
        CulcX(xPr, alf, h, x);     // X_s+1

		COUT << "Step: " << count << ENDL;
        count++;
    }

	COUT << "Num of Steps: " << count << ENDL;

	delete r;
	delete rPr;
	delete h;
	delete hPr;
	delete xPr;
}

double VectorDifSumVal(double * v1, double * v2, int n, int ProcNum)
{
    std::vector<double> resArr(omp_get_num_threads());

#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
		resArr[omp_get_thread_num()] += abs(v1[i] - v2[i]);
    }

	double res = 0;
	for each (double var in resArr)
	{
		res += var;
	}

    return res;
}

void VectorSum(double * v1, double * v2, int n, double * vRes)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        vRes[i] = v1[i] + v2[i];
    }
}

void VectorMultConst(double val, double * v1, int n, double * vRes)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        vRes[i] = v1[i] * val;
    }
}

double VectorScalarMult(double * v1, double * v2, int n, int ProcNum)
{
	std::vector<double> resArr(omp_get_num_threads());

#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
		resArr[omp_get_thread_num()] += v1[i] * v2[i];
    }

	double res = 0;
	for each (double var in resArr)
	{
		res += var;
	}

    return res;
}

void VectorMultMatrix(CRSMatrix & A, double * v1, int n, double * vRes)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        vRes[i] = 0;
        for (int j = 0; j < n; j++)
        {
            vRes[i] += A.GetValue(i, j) * v1[j];
        }
    }
}

void Vec1CopyToVec2(double * v1, double * v2, int n)
{
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		v2[i] = v1[i];
	}
}
