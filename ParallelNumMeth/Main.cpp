#include "CMatrix.h"

#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl
#define READLN int a; std::cin >> a
#define COUT_BOOL(boolVar) if (boolVar) std::cout << "True"; else std::cout << "False"

enum PNMStatus
{
    PNMStatusOk,
    PNMStatusFailed,
    PNMStatusInitError
};

void PrintVector(std::vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		COUT << vec[i] << "	";
	}

	COUT << ENDL;
}

// Блочное разложение Холецкого
// AMatrix - симетричная положительно определенная матрица
// LMatirx - нижнетреугольная с положительными элементами на диагонали
PNMStatus CholeskyBlockDecomposition(CMatrix<double> AMatrix, CMatrix<double> &LMatirx)
{
	if (AMatrix.GetSize() != LMatirx.GetSize())
	{
		return PNMStatusInitError;
	}

    for (size_t index = 0; index < pow(AMatrix.GetSize(), 2); index++)
    {
        size_t i = index / AMatrix.GetSize();
        size_t j = index % AMatrix.GetSize();

        double cuclL = 0;

        if (i == j)
        {
            double sumL = 0;

            for (size_t k = 0; k < i; k++)
            {
                sumL += pow(LMatirx[i * AMatrix.GetSize() + k], 2);
            }

            cuclL = sqrt(AMatrix[i * AMatrix.GetSize() + i] - sumL);
        }
        else if (i > j)
        {
            double sumL = 0;

            for (size_t k = 0; k < j; k++)
            {
                sumL += LMatirx[i * AMatrix.GetSize() + k] * LMatirx[j * AMatrix.GetSize() + k];
            }

            cuclL = (double)(AMatrix[index] - sumL) / (double)(LMatirx[j * AMatrix.GetSize() + j]);
        }

        LMatirx.SetMatrixCell(index, cuclL);
    }

    return PNMStatusOk;
}

PNMStatus ReverseMotion(CMatrix<double> LMatrix, std::vector<double> bVector, std::vector<double> &xVector)
{
	if (LMatrix.GetSize() == bVector.size() == xVector.size())
	{
		return PNMStatusInitError;
	}

	std::vector<double> yVector;

	for (size_t i = 0; i < LMatrix.GetSize(); i++)
	{
		double sumL = 0;

		for (size_t j = 0; j < i; j++)
		{
			sumL += LMatrix[i * LMatrix.GetSize() + j] * yVector[j];
		}

		yVector.push_back((bVector[i] - sumL) / LMatrix[i * LMatrix.GetSize() + i]);
	}

#ifdef TEST_MOD

	COUT << ENDL << "yVector : ";
	PrintVector(yVector);
	std::vector<double> refYVector = { 7.*sqrt(6.)/3., sqrt(30.)/3., sqrt(5)/2., sqrt(15.)/2. };
	COUT << "refYVector : ";
	PrintVector(refYVector);

#endif // TEST_MOD

	CMatrix<double> LTMatrix(LMatrix.GetTrMatrix());

	for (int i = LTMatrix.GetSize() - 1; i >= 0; i--)
	{
		double sumL = 0;
		int index = -1;

		for (size_t j = LTMatrix.GetSize() - 1; j > i; j--)
		{
			index = i * LTMatrix.GetSize() + j;
			sumL += LTMatrix[index] * xVector[j];
		}

		index = i * LTMatrix.GetSize() + i;

		xVector[i] = ((yVector[i] - sumL) / LTMatrix[index]);
	}

	return PNMStatusOk;
}

int main()
{
    std::vector<double> initVec = { 6,  2,   5,  1,
                                    2,  4,  -2,  4,
                                    5, -2,   9, -3,
	                                1,  4,  -3,  8  };

    CMatrix<double> AMatrix(initVec);
	COUT << "AMatrix: " << ENDL;
	AMatrix.PrintMatrix();
	COUT << "IsPositivelyDefined: ";  COUT_BOOL(AMatrix.IsPositivelyDefined()); COUT << ENDL << ENDL;

    CMatrix<double> LMatirx(AMatrix.GetSize());

    if (CholeskyBlockDecomposition(AMatrix, LMatirx) != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL;

    COUT << "LMatrix: " << ENDL;
    LMatirx.PrintMatrix();
#ifdef TEST_MOD

    COUT << ENDL << "LMatirx * LMatirx.GetTrMatrix(): " << ENDL;
	CMatrix<double> checkMatrix(LMatirx * LMatirx.GetTrMatrix());
	checkMatrix.PrintMatrix();

#endif // TEST_MOD

	/*CMatrix<double> refLMatrix({       sqrt(6.),                  0,              0,            0,
		                            sqrt(6.)/3.,       sqrt(30.)/3.,              0,            0,
		                         5.*sqrt(6.)/6., -11.*sqrt(30.)/30., 2.*sqrt(5.)/5.,            0,
	                                sqrt(6.)/6.,  11.*sqrt(30.)/30.,   sqrt(5.)/10., sqrt(15.)/2.  });
	COUT << "refLMatrix: " << ENDL;
	refLMatrix.PrintMatrix();*/

	std::vector<double> bVector = { 29, 20, 16, 32 };
	std::vector<double> xVector = {  0,  0,  0,  0 };

	COUT << ENDL << "bVector : ";
	PrintVector(bVector);
	if (ReverseMotion(LMatirx, bVector, xVector) != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL;
	COUT << "xVector : ";
	PrintVector(xVector);

	READLN;

    return 0;
}
