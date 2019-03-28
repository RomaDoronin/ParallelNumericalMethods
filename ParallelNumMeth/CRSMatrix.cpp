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
    for (int count = rowPtr[i]; count < rowPtr[i + 1]; count++)
    {
        if (colIndex[count] == j)
            return val[count];
    }

    return 0;
}

void CRSMatrix::SetValue(int i, int j, double value)
{
    if (value != 0 && GetValue(i, j) == 0)
        nz++;

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

    InsertVecPos(val, index, value);
    InsertVecPos(colIndex, index, j);

    for (int rowCount = i + 1; rowCount < n + 1; rowCount++)
        rowPtr[rowCount]++;
}

int CRSMatrix::GetN()
{
    return n;
}
