#include "CMatrix.h"

#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl
#define READLN int a; std::cin >> a

enum PNMStatus
{
    PNMStatusOk,
    PNMStatusFailed,
    PNMStatusInitError
};

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

	CMatrix<double> LTMatrix(LMatrix.GetTrMatrix());

	for (int i = LTMatrix.GetSize() - 1; i >= 0; i--)
	{
		double sumL = 0;

		for (size_t j = LTMatrix.GetSize() - 1; j > i; j--)
		{
			sumL += LMatrix[i * LMatrix.GetSize() + j] * yVector[j];
		}

		xVector[i] = ((yVector[i] - sumL) / LMatrix[i * LMatrix.GetSize() + i]);
	}

	return PNMStatusOk;
}

void PrintVector(std::vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		COUT << vec[i] << "	";
	}

	COUT << ENDL;
}

int main()
{
    std::vector<double> initVec = { 10, -3,  2,
                                    -3,  3, -2,
                                     2, -2,  7 };
	std::vector<double> bVector = { 3, 4, -14 };
	std::vector<double> xVector = {  0,  0,  0 };

    CMatrix<double> AMatrix(initVec);
    CMatrix<double> LMatirx(AMatrix.GetSize());

    if (CholeskyBlockDecomposition(AMatrix, LMatirx) != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL;

    COUT << "LMatrix: " << ENDL;
    LMatirx.PrintMatrix();
	COUT << "AMatrix: " << ENDL;
    AMatrix.PrintMatrix();
    COUT << "LMatirx * LMatirx.GetTrMatrix(): " << ENDL;
	CMatrix<double> checkMatrix(LMatirx * LMatirx.GetTrMatrix());
	checkMatrix.PrintMatrix();
    
	COUT << "bVector : ";
	PrintVector(bVector);
	if (ReverseMotion(LMatirx, bVector, xVector) != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL;
	COUT << "xVector : ";
	PrintVector(xVector);

	READLN;

    return 0;
}
