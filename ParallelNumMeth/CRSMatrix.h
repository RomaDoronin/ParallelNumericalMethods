#pragma once

#define CRS_MATRIX

#include <vector>
#include <iostream>

#define RAND(v_min, v_max) (rand() % (v_max - v_min + 1) + v_min)
#define COUT std::cout
#define CIN std::cin
#define ENDL std::endl

class CRSMatrix
{
private:
    int n; // ����������� �������
    int nz; // ����� ��������� ��������� � ����������� ������������ �������, ������� �� ���� ������� ��������� 
    std::vector<double> val; // ������ �������� ������� �� ������� 
    std::vector<int> colIndex; // ������ ������� �������� 
    std::vector<int> rowPtr; // ������ �������� ������ �����

    void InsertVecPos(std::vector<double> &vec, int index, double value);
    void InsertVecPos(std::vector<int> &vec, int index, int value);

public:
    CRSMatrix(int _n);

    double GetValue(int i, int j);

    void SetValue(int i, int j, double value);

    int GetN();
};

void GenVecWithoutNull(std::vector<double> &vec, int n, int var)
{
    for (int i = 0; i < n; i++)
    {
        int randVal = RAND(1, var);
        vec.push_back(pow(-1, RAND(0, 1)) * randVal);
    }
}

// nz - n = ������
void InitCRSMatrix(CRSMatrix &matrix, int n, int nz)
{
    std::vector<double> initVec;
    const int varNum = 9;

    // 1. ������������� ������ ����������� ((nz - n)/2) ��� 0
    GenVecWithoutNull(initVec, (nz - n) / 2, varNum);

    // 2. �������� ���������� ��� � ������� ������������
    for (int i = 0; i < initVec.size(); i++)
    {
        int indexI, indexJ;
        do
        {
            indexI = rand() % (n - 1); // RAND(0, (n - 2));
            indexJ = rand() % (n - indexI - 1) + indexI + 1; // RAND(indexI + 1, (n - 1));
        } while (matrix.GetValue(indexI, indexJ) != 0);

        matrix.SetValue(indexI, indexJ, initVec[i]);
        matrix.SetValue(indexJ, indexI, initVec[i]);
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            COUT << matrix.GetValue(i, j) << "	";
        }
        COUT << ENDL;
    }

    // 3. ��������� ��������� ������� ���������
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
}
