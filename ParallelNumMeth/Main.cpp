#include "CMatrix.h"

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
	for (size_t index = 0; index < AMatrix.GetSize() * AMatrix.GetSize(); index++)
	{
		size_t i = index / AMatrix.GetSize();
		size_t j = index % AMatrix.GetSize();

		double cuclL = 0;

		if (i == j)
		{
			double sumL = 0;

			for (size_t k = 1; k < i; k++)
			{
				sumL += pow(LMatirx[i * AMatrix.GetSize() + k], 2);
			}

			cuclL = sqrt(AMatrix[i * AMatrix.GetSize() + i] - sumL);
		}
		else if (i > j)
		{
			double sumL = 0;

			for (size_t k = 1; k < j; k++)
			{
				sumL += LMatirx[i * AMatrix.GetSize() + k] * LMatirx[j * AMatrix.GetSize() + k];
			}

			cuclL = (double)(AMatrix[index] - sumL) / (double)(LMatirx[j * AMatrix.GetSize() + j]);
		}

		LMatirx.SetMatrixCell(index, cuclL);
	}

	return PNMStatusOk;
}



int main()
{
	std::vector<double> initVec = { 2,4,9, 4,8,1, 9,1,6 };
	CMatrix<double> AMatrix(initVec);
	CMatrix<double> LMatirx(AMatrix.GetSize());

	if (CholeskyBlockDecomposition(AMatrix, LMatirx) == PNMStatusOk)
	{
		std::cout << "LMatrix: " << std::endl;
		LMatirx.PrintMatrix();

		std::cout << "AMatrix: " << std::endl;
		AMatrix.PrintMatrix();
		std::cout << "LMatirx * LMatirx.GetTrMatrix(): " << std::endl;
		(LMatirx * LMatirx.GetTrMatrix()).PrintMatrix();
		if (AMatrix == LMatirx * LMatirx.GetTrMatrix())
			std::cout << "[ OK ]";
		else
			std::cout << "[ FAILED ]";
	}
	else
	{
		std::cout << "Something goes wrong :(";
	}

	int a;
	std::cin >> a;

	return 0;
}