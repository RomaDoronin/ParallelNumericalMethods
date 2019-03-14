#include "CMatrix.h"
#include <ctime>
#include <omp.h>
#include <list>

#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl
#define READLN int a; std::cin >> a
#define COUT_BOOL(boolVar) if (boolVar) std::cout << "True"; else std::cout << "False"
#define RAND(v_min, v_max) (rand() % (v_max - v_min + 1) + v_min)
#define CHECK_STATUS(func) if (func != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL
#define CHECK_TIME(name) COUT << name << " - OK" << ENDL

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

PNMStatus CholeskyBlockDecompositionParallel(CMatrix<double> AMatrix, CMatrix<double> &LMatirx)
{
	if (AMatrix.GetSize() != LMatirx.GetSize())
	{
		return PNMStatusInitError;
	}

	// Вычисляем первый столбец

	// Параллельно вычисляем все строки, Хз пока

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

PNMStatus CholeskyBlockDecompositionWithoutCollect(CMatrix<double> AMatrix, CMatrix<double> &LMatirx)
{
	if (AMatrix.GetSize() != LMatirx.GetSize())
	{
		return PNMStatusInitError;
	}

	int size = AMatrix.GetSize();

	size_t ProcNum = 3;
	omp_set_num_threads(ProcNum);

//#pragma omp parallel for // FAILED
	for (int i = 0; i < size; i++)
	{
		AMatrix.SetMatrixCell(i * size + i, sqrt(AMatrix[i * size + i]));
#pragma omp parallel for // OK
		for (int j = i + 1; j < size; j++)
		{
			AMatrix.SetMatrixCell(j * size + i, AMatrix[j * size + i] / AMatrix[i * size + i]);
		}

#pragma omp parallel for // OK
		for (int k = i + 1; k < size; k++)
		{
			for (int j = k; j < size; j++)
			{
				AMatrix.SetMatrixCell(j * size + k, AMatrix[j * size + k] - AMatrix[j * size + i] * AMatrix[k * size + i]);
			}
		}
	}


	for (int i = 0; i < pow(size, 2); i++)
	{
		if (i % size <= i / size)
			LMatirx.SetMatrixCell(i, AMatrix[i]);
		else
			LMatirx.SetMatrixCell(i, 0);
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

std::vector<double> GenVec(int size, int var)
{
	std::vector<double> res;

	for (int i = 0; i < size; i++)
	{
		res.push_back(RAND(-var, var));
	}

	return res;
}

// Positively Defined Simetric Matrix
// koef - The power ratio of the main diagonal
PNMStatus GeneratePDSM(int size, double koef, int var, std::vector<double> &initVec)
{
	std::vector<double> vec(pow(size, 2));

	for (int i = 0; i < pow(size, 2); i++)
	{
		int rowIndex = i / size;
		int colIndex = i % size;
		
		if (rowIndex < colIndex)
		{
			vec[rowIndex * size + colIndex] = RAND(-var, var);
			vec[colIndex * size + rowIndex] = vec[rowIndex * size + colIndex];
		}
		else if (rowIndex == colIndex)
		{
			vec[i] = RAND(-var, var);
		}
	}

	int sum = 0;

	for (int i = 0; i < pow(size, 2); i++)
	{
		int rowIndex = i / size;
		int colIndex = i % size;

		if (rowIndex != colIndex)
			sum += abs(vec[i]);

		if (colIndex == size - 1)
		{
			while (vec[rowIndex * size + rowIndex] < sum * koef)
			{
				vec[rowIndex * size + rowIndex] += RAND(1, var);
			}

			sum = 0;
		}
	}

	initVec = vec;

	return PNMStatusOk;
}

PNMStatus PRKK(std::vector<double> res1, std::vector<double> res2, double accuracy)
{
	if (res1.size() != res2.size())
		return PNMStatusFailed;

	for (int i = 0; i < res1.size(); i++)
	{
		if (abs(res1[i] - res2[i]) > accuracy)
			return PNMStatusFailed;
	}

	return PNMStatusOk;
}

std::vector<double> GetCopyVector(std::vector<double> vec)
{
	std::vector<double> c_vec;

	for each (double var in vec)
	{
		c_vec.push_back(var);
	}

	return vec;
}

int main()
{
	srand(time(0));
	int size = 0;

	/*std::list<std::vector<double>> dataList;

	for (size = 1000; size < 5001;)
	{
		std::vector<double> initVec(pow(size, 2));
		CHECK_STATUS(GeneratePDSM(size, 1, 9, initVec));
		CMatrix<double> AMatrix(initVec);
		if (AMatrix.IsPositivelyDefined())
		{
			COUT << "PUSH" << ENDL;
			dataList.push_back(initVec);
			size += 1000;
		}
	}*/

	for (int co = 0; co < 3; co++)
	{
		for (size = 1000; size < 5001; size += 1000)
		{
			double start_time = omp_get_wtime();
			double curr_time = start_time;

			std::vector<double> initVec(pow(size, 2));
			CHECK_STATUS(GeneratePDSM(size, 1, 9, initVec));
			CHECK_TIME("Gen");
			CMatrix<double> AMatrix(initVec /*{ 2,-1,1, -1,4,2, 1,2,6 }*/);

			//COUT << "AMatrix: " << ENDL;
			//AMatrix.PrintMatrix();
			//COUT << "IsPositivelyDefined: ";  COUT_BOOL(AMatrix.IsPositivelyDefined()); COUT << ENDL << ENDL;

			std::vector<double> ref_xVector = GenVec(size, 20);
			//COUT << "Ref xVector : ";
			//PrintVector(ref_xVector);
			std::vector<double> bVector = AMatrix * ref_xVector;
			//COUT << "bVector : ";
			//PrintVector(bVector);


			CMatrix<double> LMatirx(AMatrix.GetSize());

			//CHECK_STATUS(CholeskyBlockDecomposition(AMatrix, LMatirx));
			CHECK_STATUS(CholeskyBlockDecompositionWithoutCollect(AMatrix, LMatirx));
			CHECK_TIME("Decomposition");

			//COUT << ENDL << "LMatrix: " << ENDL;
			//LMatirx.PrintMatrix();

#ifdef TEST_MOD

			COUT << ENDL << "LMatirx * LMatirx.GetTrMatrix(): " << ENDL;
			CMatrix<double> checkMatrix(LMatirx * LMatirx.GetTrMatrix());
			checkMatrix.PrintMatrix();

#endif // TEST_MOD

			std::vector<double> xVector(size);
			CHECK_STATUS(ReverseMotion(LMatirx, bVector, xVector));
			CHECK_TIME("Reverse");

			//COUT << ENDL << "xVector : ";
			//PrintVector(xVector);

			CHECK_STATUS(PRKK(ref_xVector, xVector, 0.00001));
			CHECK_TIME("PRKK");

			COUT << "====================" << ENDL;
			COUT << "Size: " << size << " | Time: " << omp_get_wtime() - start_time << ENDL;
			COUT << "====================" << ENDL;
		}

		COUT << ".............................................................." << ENDL;
	}

	READLN;

    return 0;
}
