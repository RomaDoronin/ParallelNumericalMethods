#include "CMatrix.h"
#include <ctime>
#include <omp.h>
#include <list>

#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl
#define READLN std::cout << "Input any key..." << ENDL; int a; std::cin >> a
#define COUT_BOOL(boolVar) if (boolVar) std::cout << "True"; else std::cout << "False"
#define RAND(v_min, v_max) (rand() % (v_max - v_min + 1) + v_min)
#define CHECK_STATUS(func) if (func != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL
#define CHECK_TIME(name) COUT << name << " - OK" << ENDL

#define ACCURACY 0.00001

//#define CHOLESKY_DECOMPOSITION
#define CONJUGATE_GRADIENT_METHOD

//////////////////////////////////////////////////
///////////////// Cholesky decomposition
//////////////////////////////////////////////////
// Реализовать блочное разложение Холецкого для симметричной положительно определенной матрицы А 
// AMatrix - симетричная положительно определенная матрица
// LMatirx - нижнетреугольная с положительными элементами на диагонали
void CholeskyBlockDecompositionWithoutCollect(CMatrix<double> AMatrix, CMatrix<double> &LMatirx, int ProcNum)
{
	int size = AMatrix.GetSize();

	omp_set_num_threads(ProcNum);


	for (int i = 0; i < size; i++)
	{
		// Извлекаем корень из i-ого элемента на главной диагонали
		// A[i][i] = sqtr(A[i][i])
		AMatrix.SetMatrixCell(i * size + i, sqrt(AMatrix[i * size + i]));

#ifdef FILE_WORK
		AMatrix.WriteMatrixInFile("tread_0.txt", "A[" + std::to_string(i) + "][" + std::to_string(i) + "] = sqtr(A[" + std::to_string(i) + "][" + std::to_string(i) + "])"
			+ " i = " + std::to_string(i));
#endif

		// Делим все элеметны i-ого столбца ниже главной диагонали на i-ый элемент главной диагонали 
		// A[j][i] /= A[i][i]
#pragma omp parallel for // OK
		for (int j = i + 1; j < size; j++)
		{
			AMatrix.SetMatrixCell(j * size + i, AMatrix[j * size + i] / AMatrix[i * size + i]);
			
#ifdef FILE_WORK
			std::string filename = "tread_" + std::to_string(omp_get_thread_num()) + ".txt";
			AMatrix.WriteMatrixInFile(filename, "A[" + std::to_string(j) + "][" + std::to_string(i) + "] = A[" + std::to_string(j) + "][" + std::to_string(i) + "] / A[" + std::to_string(i) + "][" + std::to_string(i) + "]"
				+ " i = " + std::to_string(i) + " j = " + std::to_string(j));
#endif
		}

		// A[j][k] = A[j][k] - A[j][i] * A[k][i]
//#pragma omp parallel for // OK
		for (int k = i + 1; k < size; k++)
		{
#pragma omp parallel for // OK
			for (int j = k; j < size; j++)
			{
				AMatrix.SetMatrixCell(j * size + k, AMatrix[j * size + k] - AMatrix[j * size + i] * AMatrix[k * size + i]);

#ifdef FILE_WORK
				std::string filename = "tread_" + std::to_string(omp_get_thread_num()) + ".txt";
				AMatrix.WriteMatrixInFile(filename, "A[" + std::to_string(j) + "][" + std::to_string(k) + "] = A[" + std::to_string(j) + "][" + std::to_string(k) + "] - A[" + std::to_string(j) + "][" + std::to_string(i) + "] * A[" + std::to_string(k) + "][" + std::to_string(i) + "]"
					+ " i = " + std::to_string(i) + " j = " + std::to_string(j) + " k = " + std::to_string(k));
#endif
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
}

void ReverseMotion(CMatrix<double> LMatrix, std::vector<double> bVector, std::vector<double> &xVector)
{
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
	std::vector<double> refYVector = { 7.*sqrt(6.) / 3., sqrt(30.) / 3., sqrt(5) / 2., sqrt(15.) / 2. };
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
}


//////////////////////////////////////////////////
///////////////// Conjugate gradient method
//////////////////////////////////////////////////
// Реализовать метод сопряженных градиентов для решения СЛАУ Ax=b с симметричной положительно определенной разреженной матрицей A и плотным вектором b. 
//
//

//////////////////////////////////////////////////
///////////////// Softgrader
//////////////////////////////////////////////////
void Cholesky_Decomposition(double * A, double * L, int n)
{
    for (int i = 0; i < n; i++)
    {
        // Извлекаем корень из i-ого элемента на главной диагонали
        // A[i][i] = sqtr(A[i][i])
        A[i * n + i] = sqrt(A[i * n + i]);

        // Делим все элеметны i-ого столбца ниже главной диагонали на i-ый элемент главной диагонали 
        // A[j][i] /= A[i][i]
#pragma omp parallel for
        for (int j = i + 1; j < n; j++)
        {
            A[j * n + i] /= A[i * n + i];
        }

        // A[j][k] = A[j][k] - A[j][i] * A[k][i]
        for (int k = i + 1; k < n; k++)
        {
#pragma omp parallel for // OK
            for (int j = k; j < n; j++)
            {
                A[j * n + k] -= A[j * n + i] * A[k * n + i];
			}
        }
    }

	for (int i = 0; i < pow(n, 2); i++)
	{
		if (i % n <= i / n)
			L[i] = A[i];
		else
			L[i] = 0;
	}
}


//////////////////////////////////////////////////
///////////////// Local function
//////////////////////////////////////////////////
void PrintVector(std::vector<double> vec)
{
	for (int i = 0; i < vec.size(); i++)
	{
		COUT << vec[i] << "	";
	}

	COUT << ENDL;
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
void GeneratePDSM(int size, double koef, int var, std::vector<double> &initVec)
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
}

bool PRKK(std::vector<double> res1, std::vector<double> res2, double accuracy)
{
	if (res1.size() != res2.size())
		return false;

	for (int i = 0; i < res1.size(); i++)
	{
		if (abs(res1[i] - res2[i]) > accuracy)
			return false;
	}

	return true;
}

std::vector<double> GetCopyVector(std::vector<double> vec)
{
	std::vector<double> c_vec;

	for (const auto &var: vec)
	{
		c_vec.push_back(var);
	}

	return vec;
}

bool CompareMatrixVsArray(CMatrix<double> matrix, double *arr)
{
	for (int i = 0; i < matrix.GetSize(); i++)
	{
		if (matrix[i] != arr[i])
		{
			return false;
		}
	}

	return true;
}


//////////////////////////////////////////////////
///////////////// Main function
//////////////////////////////////////////////////
int main()
{
#ifdef CHOLESKY_DECOMPOSITION

#ifdef TESTPROG
	srand(time(0));

	int ProcNum = 4;

#ifdef FILE_WORK
	for (int i = 0; i < ProcNum; i++)
	{
		std::ofstream fout("tread_" + std::to_string(i) + ".txt", std::ios_base::trunc);
		fout.close();
	}
#endif

	for (int size = 354; size <= 2000; size += 139)
	{
		COUT << "=================================== SIZE : " << size << ENDL;

		std::vector<double> initVec(pow(size, 2));
		CHECK_STATUS(GeneratePDSM(size, 1, 9, initVec));
		CMatrix<double> AMatrixRef(initVec);
		CHECK_TIME("Gen");

		for (int countTreadNum = 1; countTreadNum < ProcNum + 1; countTreadNum++)
		{
			COUT << "=================================== COUNT TREAD NUM : " << countTreadNum << ENDL;

			for (int count = 0; count < 2; count++)
			{
				COUT << "=================================== COUNT : " << count << ENDL;
				
				CMatrix<double> AMatrix(AMatrixRef);
				std::vector<double> ref_xVector = GenVec(AMatrix.GetSize(), 20);
				std::vector<double> bVector = AMatrix * ref_xVector;
				CMatrix<double> LMatirx(AMatrix.GetSize());

				double start_time = omp_get_wtime();

				CHECK_STATUS(CholeskyBlockDecompositionWithoutCollect(AMatrix, LMatirx, countTreadNum));
				CHECK_TIME("Decomposition");

				COUT << "====================" << ENDL;
				COUT << "Time: " << omp_get_wtime() - start_time << ENDL;
				COUT << "====================" << ENDL;

				std::vector<double> xVector(size);
				CHECK_STATUS(ReverseMotion(LMatirx, bVector, xVector));
				CHECK_TIME("Reverse");

				CHECK_STATUS(PRKK(ref_xVector, xVector, ACCURACY));
				CHECK_TIME("PRKK");
			}
		}
	}
#endif

	int size = 3;
	int ProcNum = 1;

	std::vector<double> initVec(pow(size, 2));
	CHECK_STATUS(GeneratePDSM(size, 1, 9, initVec));

	CMatrix<double> AMatrix(initVec);
	std::vector<double> ref_xVector = GenVec(AMatrix.GetSize(), 20);
	std::vector<double> bVector = AMatrix * ref_xVector;
	CMatrix<double> LMatirx(AMatrix.GetSize());

	CHECK_STATUS(CholeskyBlockDecompositionWithoutCollect(AMatrix, LMatirx, ProcNum));

	// -------------------------------------------------------------------------------
	int n = size;
	double *A = new double[pow(n, 2)];
	for (int i = 0; i < pow(n,2); i++)
	{
		A[i] = initVec[i];
	}
	double *L = new double[pow(n, 2)];

	omp_set_num_threads(ProcNum);
	Cholesky_Decomposition(A, L, n);

	if (CompareMatrixVsArray(LMatirx, L))
		COUT << "OK";
	else
		COUT << "FAILED";
#endif

#ifdef CONJUGATE_GRADIENT_METHOD
	


#endif

	READLN;

    return 0;
}
