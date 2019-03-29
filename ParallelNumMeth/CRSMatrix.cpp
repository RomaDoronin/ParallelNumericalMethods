#include "CRSMatrix.h"

void CRSMatrix::InsertVecPos(std::vector<double> &vec, int index, double value)
{
    std::vector<double> vecRes;

    for (int i = 0; i < vec.size(); i++)
    {
        if (i == index)
            vecRes.push_back(value);

        vecRes.push_back(vec[i]);
    }

    if (vec.size() <= index)
    {
        vecRes.push_back(value);
    }

    vec = vecRes;
}

void CRSMatrix::InsertVecPos(std::vector<int> &vec, int index, int value)
{
    std::vector<int> vecRes;

    for (int i = 0; i < vec.size(); i++)
    {
        if (i == index)
            vecRes.push_back(value);

        vecRes.push_back(vec[i]);
    }

    if (vec.size() <= index)
    {
        vecRes.push_back(value);
    }

    vec = vecRes;
}

CRSMatrix::CRSMatrix(int _n)
{
    n = _n;
    nz = 0;
    for (int i = 0; i < n; i++)
    {
        rowPtr.push_back(0);
    }

    rowPtr.push_back(0);
}

double CRSMatrix::GetValue(int i, int j)
{
	if (j < i)
	{
		int tmp = i;
		i = j;
		j = tmp;
	}

    for (int count = rowPtr[i]; count < rowPtr[i + 1]; count++)
    {
        if (colIndex[count] == j)
            return val[count];
    }

    return 0;
}

void CRSMatrix::SetValue(int i, int j, double value)
{
    /*if (value != 0 && GetValue(i, j) == 0)
        nz++;*/

    int index = rowPtr[i + 1];

    for (int count = rowPtr[i]; count < rowPtr[i + 1]; count++)
    {
        if (j < colIndex[count])
        {
            index = count;
            break;
        }
        else if (j == colIndex[count])
        {
            val[count] = value;
            colIndex[count] = j;
            return;
        }
    }

	val.insert(val.begin() + index, value);
    //InsertVecPos(val, index, value);
	colIndex.insert(colIndex.begin() + index, j);
    //InsertVecPos(colIndex, index, j);

    for (int rowCount = i + 1; rowCount < n + 1; rowCount++)
        rowPtr[rowCount]++;
}

int CRSMatrix::GetN()
{
    return n;
}

void GenVecWithoutNull(std::vector<double> &vec, int n, int var)
{
//#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        int randVal = RAND(1, var);
        vec.push_back(pow(-1, RAND(0, 1)) * randVal);
    }
}

void InitCRSMatrix(CRSMatrix &matrix, int n, int nz)
{
    std::vector<double> initVec;
    const int varNum = 9;

	COUT << "NZ: " << nz << ENDL;

    // 1. Сгенерировать вектор размерности ((nz - n)/2) без 0
    GenVecWithoutNull(initVec, (nz - n) / 2, varNum);

    // 2. Рандомно разместить его в верхнем треугольнике
//#pragma omp parallel for
    for (int i = 0; i < initVec.size(); i++)
    {
		//COUT << "i: " << i << ENDL;

        int indexI, indexJ;
        //do
        //{
            indexI = rand() % (n - 1); // RAND(0, (n - 2));
            indexJ = rand() % (n - indexI - 1) + indexI + 1; // RAND(indexI + 1, (n - 1));
        //} while (matrix.GetValue(indexI, indexJ) != 0);

        matrix.SetValue(indexI, indexJ, initVec[i]);
    }

    // 3. Заполнить Правильно главную диагональ
//#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < n; j++)
        {
            sum += abs(matrix.GetValue(i, j));
        }

        if (sum == 0)
            sum = RAND(1, varNum);

        matrix.SetValue(i, i, sum);
    }

	if (n < 5)
	{
		PrintCRSMatrix(matrix);
	}
}

void PrintCRSMatrix(CRSMatrix &matrix)
{
	COUT << ENDL << "CRSMatrix: " << ENDL;
	for (int i = 0; i < matrix.GetN(); i++)
	{
		for (int j = 0; j < matrix.GetN(); j++)
		{
			COUT << matrix.GetValue(i, j) << "	";
		}
		COUT << ENDL;
	}
}
