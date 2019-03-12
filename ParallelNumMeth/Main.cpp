#include "CMatrix.h"

#define COUT std::cout
#define ENDL std::endl

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


int main()
{
    std::vector<double> initVec = { 10, -3,  2,
                                    -3,  3, -2,
                                     2, -2,  7 };
    CMatrix<double> AMatrix(initVec);
    CMatrix<double> LMatirx(AMatrix.GetSize());

    if (CholeskyBlockDecomposition(AMatrix, LMatirx) == PNMStatusOk)
    {
        std::cout << "LMatrix: " << std::endl;
        LMatirx.PrintMatrix();

        std::cout << "AMatrix: " << std::endl;
        AMatrix.PrintMatrix();
        std::cout << "LMatirx * LMatirx.GetTrMatrix(): " << std::endl;
		CMatrix<double> checkMatrix(LMatirx * LMatirx.GetTrMatrix());
		checkMatrix.PrintMatrix();
    }
    else
    {
        std::cout << "Something goes wrong :(";
    }

    int a;
    std::cin >> a;

    return 0;
}
