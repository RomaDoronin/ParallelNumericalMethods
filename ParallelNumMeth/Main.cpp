#ifndef _CMATRIX_
#include "CMatrix.h"
#endif
#ifndef _CONJUGATEGRADIENTMETHOD_
#include "ConjugateGradientMethod.h"
#endif
#ifndef _CRSMATRIX_
#include "CRSMatrix.h"
#endif
#ifndef _CTIME_
#include <ctime>
#endif
#ifndef _LIST_
#include <list>
#endif
#ifndef _FSTREAM_
#include <fstream>
#endif

#ifndef COUT
#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl
#endif // COUT

#define READLN std::cout << "Input any key..." << ENDL; int a; std::cin >> a
#define COUT_BOOL(boolVar) if (boolVar) std::cout << "True"; else std::cout << "False"
#ifndef RAND
#define RAND(v_min, v_max) (rand() % (v_max - v_min + 1) + v_min)
#endif
#define CHECK_STATUS(func) if (func != PNMStatusOk) COUT << "Something goes wrong :(" << ENDL
#define CHECK_TIME(name) COUT << name << " - OK" << ENDL

#define CHOLESKY_DECOMPOSITION
//#define CONJUGATE_GRADIENT_METHOD

/**/
#define ELEM_NUM 10000
#define PROC_NUM 8
/**/

void PrintVector(double * vec, int rowSize, int colSize);

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

        // Делим все элеметны i-ого столбца ниже главной диагонали на i-ый элемент главной диагонали 
        // A[j][i] /= A[i][i]
#pragma omp parallel for
        for (int j = i + 1; j < size; j++)
        {
            AMatrix.SetMatrixCell(j * size + i, AMatrix[j * size + i] / AMatrix[i * size + i]);
        }

        // A[j][k] = A[j][k] - A[j][i] * A[k][i]
        for (int k = i + 1; k < size; k++)
        {
#pragma omp parallel for
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
///////////////// Softgrader
//////////////////////////////////////////////////
// r - Место раздела матрицы
double * GetSubMatrix(double * A, int n, int row1, int row2, int col1, int col2)
{
    double * subMatrix = new double[(row2 - row1) * (col2 - col1)];
    int count = 0;

    for (int i = row1; i < row2; i++)
    {
        for (int j = col1; j < col2; j++)
        {
            subMatrix[count] = A[i * n + j];
            count++;
        }
    }

    return subMatrix;
}

// A22 - L21*L21T
double * WaveOperation(double * A22, double * L21, int n, int r)
{
    int size = n - r;
    double * res = new double[size * size];

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            int multVec = 0;
            for (int count = 0; count < r; count++)
            {
                //std::cout << "L21[i * r + count] : " << L21[i * r + count] << " | L21[j * r + count] : " << L21[j * r + count] << std::endl;
                multVec += L21[i * r + count] * L21[j * r + count];
            }
            //std::cout << "A22[i * size + j] : " << A22[i * size + j] << " | multVec : " << multVec << std::endl;
            res[i * size + j] = A22[i * size + j] - multVec;
        }
    }

    return res;
}

// Решается множество систем линейных уравнений: L21 * L11T = A21 -> L11 * L21T = A21T
void CulcL21(double * A21, double * L21, double * L11, int n, int r)
{
    int size = n - r;

    // Преведение всей системы к треугольному виду
    /*for (int mainDiagCount = 0; mainDiagCount < r; mainDiagCount++)
    {
    for (int rowCount = mainDiagCount + 1; rowCount < r; rowCount++)
    {
    int d = L11[rowCount * r + mainDiagCount] / L11[mainDiagCount * r + mainDiagCount];

    for (int colCount = mainDiagCount; colCount < r; colCount++)
    {
    L11[rowCount * r + colCount] -= d * L11[mainDiagCount * r + colCount];
    }

    for (int count = 0; count < size; count++)
    {
    //A21[rowCount * size + count] -= d * A21[mainDiagCount * r + count];
    // A21T
    A21[count * size + rowCount] -= d * A21[count * r + mainDiagCount];
    }
    }
    }*/

    // Обратный ход для каждой правой части
    for (int numRightPart = 0; numRightPart < size; numRightPart++)
    {
        for (int colCount = 0; colCount < r; colCount++)
        {
            int sum = 0;

            for (int rowCount = 0; rowCount < colCount; rowCount++)
            {
                // L11T
                sum += L11[colCount * r + rowCount] * L21[numRightPart * r + rowCount];
            }

            //std::cout << "A21[numRightPart * r + colCount] : " << A21[numRightPart * r + colCount] << " | sum : " << sum << " | L11[colCount * r + colCount] : " << L11[colCount * r + colCount] << std::endl;
            L21[numRightPart * r + colCount] = (A21[numRightPart * r + colCount] - sum) / L11[colCount * r + colCount];
        }
    }
}

void MergeLMatrix(double * L, double * L11, double * L21, double * L22, int n, int blockSize)
{
    int count = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i < blockSize && j < blockSize)
            {
                L[count] = L11[i * blockSize + j];
            }
            else if (i < blockSize && j >= blockSize)
            {
                L[count] = 0;
            }
            else if (i >= blockSize && j < blockSize)
            {
                L[count] = L21[(i - blockSize) * blockSize + j];
            }
            else if (i >= blockSize && j >= blockSize)
            {
                L[count] = L22[(i - blockSize) * (n - blockSize) + (j - blockSize)];
            }

            count++;
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
        // Идем по оставшуся нижнему треугольнику
        for (int k = i + 1; k < n; k++)
        {
#pragma omp parallel for
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
    // Размер блока
    int blockSize = 3;

    double * A11 = GetSubMatrix(A, n, 0, blockSize, 0, blockSize); PrintVector(A11, blockSize, blockSize);
    double * A21 = GetSubMatrix(A, n, blockSize, n, 0, blockSize); PrintVector(A21, n - blockSize, blockSize);
    double * A22 = GetSubMatrix(A, n, blockSize, n, blockSize, n); PrintVector(A22, n - blockSize, n - blockSize);

    double * L11 = new double[blockSize * blockSize];
    double * L21 = new double[(n - blockSize) * blockSize];
    double * L22 = new double[(n - blockSize) * (n - blockSize)];

    Cholesky_Decomposition_Func(A11, L11, blockSize); PrintVector(L11, blockSize, blockSize);

    CulcL21(A21, L21, L11, n, blockSize);
    
    double * _A22 = WaveOperation(A22, L21, n, blockSize); PrintVector(_A22, n - blockSize, n - blockSize);
    Cholesky_Decomposition_Func(_A22, L22, n - blockSize); PrintVector(L22, n - blockSize, n - blockSize);

    MergeLMatrix(L, L11, L21, L22, n, blockSize); PrintVector(L, n, n);

    delete [] A11;
    delete [] A21;
    delete [] A22;

    delete [] L11;
    delete [] L21;
    delete [] L22;

    delete [] _A22;
}

//////////////////////////////////////////////////
///////////////// Local function
//////////////////////////////////////////////////
void PrintVector(double * vec, int rowSize, int collSize)
{
    for (int i = 0; i < rowSize; i++)
    {
        for (int j = 0; j < collSize; j++)
        {
            COUT << vec[i * collSize + j] << "	";
        }
        COUT << ENDL;
    }

    COUT << ENDL;
}

std::vector<double> GenVec(int size, int var)
{
    std::vector<double> res;

    for (int i = 0; i < size; i++)
    {
        res.push_back(0);

        do {
            res[i] = RAND(-var, var);
        } while (res[i] == 0);
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
            vec[rowIndex * size + colIndex] = RAND(0, var * 2);
            vec[colIndex * size + rowIndex] = vec[rowIndex * size + colIndex];
        }
        else if (rowIndex == colIndex)
        {
            vec[i] = RAND(0, var * 2);
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

bool PRKK(double * res1, double * res2, int n, double accuracy)
{
    for (int i = 0; i < n; i++)
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
		double val = matrix[i];
        if (val != arr[i])
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
    srand(time(0));

#ifdef TESTPROG

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

    int size = 9;
    int ProcNum = 1;

    std::vector<double> initVec(pow(size, 2));
    GeneratePDSM(size, 1, 9, initVec);

    CMatrix<double> AMatrix(initVec);
    std::vector<double> ref_xVector = GenVec(AMatrix.GetSize(), 20);
    std::vector<double> bVector = AMatrix * ref_xVector;
    CMatrix<double> LMatirx(AMatrix.GetSize());

    CholeskyBlockDecompositionWithoutCollect(AMatrix, LMatirx, ProcNum);
	LMatirx.PrintMatrix();
    // -------------------------------------------------------------------------------
    int n = size;
    double *A = new double[pow(n, 2)];
    for (int i = 0; i < pow(n,2); i++)
    {
        A[i] = initVec[i];
    }
    PrintVector(A, n, n);
    double *L = new double[pow(n, 2)];

    omp_set_num_threads(ProcNum);
    Cholesky_Decomposition(A, L, n);

    if (CompareMatrixVsArray(LMatirx, L))
        COUT << "OK";
    else
        COUT << "FAILED";

    COUT << ENDL;

    READLN;

    return 0;
}
