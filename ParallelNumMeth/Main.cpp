#include <iostream>
#include <omp.h>

#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl


#define READLN std::cout << "Input any key..." << ENDL; int a; std::cin >> a
#define COUT_BOOL(boolVar) if (boolVar) std::cout << "True"; else std::cout << "False"
#ifndef RAND
#define RAND(v_min, v_max) (rand() % (v_max - v_min + 1) + v_min)
#endif
#define CHECK_STATUS(func) if (func != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL
#define CHECK_TIME(name) COUT << name << " - OK" << ENDL

//////////////////////////////////////////////////
///////////////// Cholesky decomposition
//////////////////////////////////////////////////
///////////////// Softgrader
//////////////////////////////////////////////////
// r - Место раздела матрицы
double * GetSubMatrix(double * A, int n, int row1, int row2, int col1, int col2)
{
	int sizeI = row2 - row1;
	int sizeJ = col2 - col1;

    double * subMatrix = new double[sizeI * sizeJ];

#pragma omp parallel for
    for (int i = row1; i < row2; i++)
    {
        for (int j = col1; j < col2; j++)
        {
            subMatrix[(i - row1) * sizeJ + (j - col1)] = A[i * n + j];
        }
    }

    return subMatrix;
}

// A22 - L21*L21T
double * WaveOperation(double * A22, double * L21, int n, int r)
{
    int size = n - r;
    double * res = new double[size * size];

#pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            double sum = 0;

            for (int count = 0; count < r; count++)
            {
                sum += L21[i * r + count] * L21[j * r + count];
            }

            res[i * size + j] = A22[i * size + j] - sum;
        }
    }

    return res;
}

// Решается множество систем линейных уравнений: L21 * L11T = A21 -> L11 * L21T = A21T
void CulcL21(double * A21, double * L21, double * L11, int n, int r)
{
    int size = n - r;

// Обратный ход обычный
#pragma omp parallel for
    for (int blockCount = 0; blockCount < size; blockCount++)
    {
        for (int colCount = 0; colCount < r; colCount++)
        {
            double sum = 0;

            for (int rowCount = 0; rowCount < colCount; rowCount++)
            {
                sum += L11[colCount * r + rowCount] * L21[blockCount * r + rowCount];
            }

            // L11T
            L21[blockCount * r + colCount] = (A21[blockCount * r + colCount] - sum) / L11[colCount * r + colCount];
        }
    }
}

void MergeL11(double * L, double * L11, int n, int blockSize)
{
#pragma omp parallel for
	for (int i = 0; i < blockSize; i++)
	{
		for (int j = 0; j < blockSize; j++)
		{
			L[i * n + j] = L11[i * blockSize + j];
		}
	}
}

void MergeL12(double * L, int n, int blockSize)
{
#pragma omp parallel for
	for (int i = 0; i < blockSize; i++)
	{
		for (int j = blockSize; j < n; j++)
		{
			L[i * n + j] = 0;
		}
	}
}

void MergeL21(double * L, double * L21, int n, int blockSize)
{
#pragma omp parallel for
	for (int i = blockSize; i < n; i++)
	{
		for (int j = 0; j < blockSize; j++)
		{
			L[i * n + j] = L21[(i - blockSize) * blockSize + j];
		}
	}
}

void MergeL22(double * L, double * L22, int n, int blockSize)
{
#pragma omp parallel for
	for (int i = blockSize; i < n; i++)
	{
		for (int j = blockSize; j < n; j++)
		{
			L[i * n + j] = L22[(i - blockSize) * (n - blockSize) + (j - blockSize)];
		}
	}
}

void Cholesky_Decomposition_Func(double * A, double * L, int n)
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

void Cholesky_Decomposition(double * A, double * L, int n)
{
    int blockSize = 125;

    if (n > blockSize)
    {
        double * A11 = GetSubMatrix(A, n, 0, blockSize, 0, blockSize);
        double * L11 = new double[blockSize * blockSize];
        Cholesky_Decomposition_Func(A11, L11, blockSize);
		MergeL11(L, L11, n, blockSize);
        delete[] A11;

        double * A21 = GetSubMatrix(A, n, blockSize, n, 0, blockSize);
        double * L21 = new double[(n - blockSize) * blockSize];
        CulcL21(A21, L21, L11, n, blockSize);
		MergeL21(L, L21, n, blockSize);
		delete[] L11;
        delete[] A21;

        double * A22 = GetSubMatrix(A, n, blockSize, n, blockSize, n);
        double * L22 = new double[(n - blockSize) * (n - blockSize)];
        double * _A22 = WaveOperation(A22, L21, n, blockSize);
        delete[] A22;
		delete[] L21;
        Cholesky_Decomposition(_A22, L22, n - blockSize);
        delete[] _A22;
		MergeL22(L, L22, n, blockSize);
		delete[] L22;

		MergeL12(L, n, blockSize);
    }
    else
    {
        Cholesky_Decomposition_Func(A, L, n);
    }
}

//////////////////////////////////////////////////
///////////////// Local function
//////////////////////////////////////////////////
void PrintMatrix(double * matrix, int rowSize, int collSize)
{
    COUT << "Matrix" << ENDL;

    for (int i = 0; i < rowSize; i++)
    {
        for (int j = 0; j < collSize; j++)
        {
            COUT << ((int)(matrix[i * collSize + j] * 100)) / 100. << "	";
        }
        COUT << ENDL;
    }

    COUT << ENDL;
}

double* GenVec(int size, int var)
{
    double * res = new double[size * size];

    for (int i = 0; i < size * size; i++)
    {
        res[i] = 0;

        do {
            res[i] = RAND(-var, var);
        } while (res[i] == 0);
    }

    return res;
}

// Positively Defined Simetric Matrix
// koef - The power ratio of the main diagonal
void GeneratePDSM(int size, double koef, int var, double * initMatrix)
{
    for (int i = 0; i < pow(size, 2); i++)
    {
        int rowIndex = i / size;
        int colIndex = i % size;
        
        if (rowIndex < colIndex)
        {
            initMatrix[rowIndex * size + colIndex] = RAND(0, var * 2);
            initMatrix[colIndex * size + rowIndex] = initMatrix[rowIndex * size + colIndex];
        }
        else if (rowIndex == colIndex)
        {
            initMatrix[i] = RAND(0, var * 2);
        }
    }

    int sum = 0;

    for (int i = 0; i < pow(size, 2); i++)
    {
        int rowIndex = i / size;
        int colIndex = i % size;

        if (rowIndex != colIndex)
            sum += abs(initMatrix[i]);

        if (colIndex == size - 1)
        {
            while (initMatrix[rowIndex * size + rowIndex] < sum * koef)
            {
                initMatrix[rowIndex * size + rowIndex] += RAND(1, var);
            }

            sum = 0;
        }
    }
}

bool PRKK(double * res1, double * res2, int size, double accuracy)
{
    for (int i = 0; i < size * size; i++)
    {
        if (abs(res1[i] - res2[i]) > accuracy)
            return false;
    }

    return true;
}

void CopyMatrix(double * matrix, double * copy_matrix, int size)
{
    for (int i = 0; i < size * size; i++)
    {
        copy_matrix[i] = matrix[i];
    }
}

//double* MultSquareMatrix(double * mat1, double * mat2, int size)
//{
//	double * res = new double[size*size];
//
//	for (int i = 0; i < size; i++)
//	{
//		for (int j = 0; j < size; j++)
//		{
//			res[i * size + j] = 0;
//
//			for (int v = 0; v < size; v++)
//			{
//				res[i * size + j] += mat1[i * size + v] * mat2[v * size + j];
//			}
//		}
//	}
//
//	return res;
//}

//////////////////////////////////////////////////
///////////////// Main function
//////////////////////////////////////////////////
int main()
{
    int size = 1000;
	COUT << "Matrix size: " << size << ENDL;

	double start_time = 0.0;

    int ProcNum = 5;
    omp_set_num_threads(ProcNum);
	COUT << "Process Number: " << ProcNum << ENDL;

    double * A = new double[size * size];
    GeneratePDSM(size, 1, 10, A);

    // -------------------------------------------------------------------------------
    double *LRef = new double[size * size];
    double * c_A = new double[size * size];
    CopyMatrix(A, c_A, size);
	start_time = omp_get_wtime();
    Cholesky_Decomposition_Func(c_A, LRef, size);
	COUT << "Reference Time: " << omp_get_wtime() - start_time << ENDL;

    // -------------------------------------------------------------------------------
    double *L = new double[size * size];
	start_time = omp_get_wtime();
    Cholesky_Decomposition(A, L, size);
	COUT << "Origin Time: " << omp_get_wtime() - start_time << ENDL;

    // -------------------------------------------------------------------------------
	COUT << "PRKK: ";
    if (PRKK(LRef, L, size, 0.00000001))
    {
        COUT << "OK";
    }
    else
    {
        COUT << "FAILED";
    }

    COUT << ENDL;

    READLN;

    return 0;
}
