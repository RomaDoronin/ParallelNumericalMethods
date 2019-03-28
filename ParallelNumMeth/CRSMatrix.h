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
    int n; // Размерность матрицы
    int nz; // Число ненулевых элементов в разреженной симметричной матрице, лежащих не ниже главной диагонали 
    std::vector<double> val; // Массив значений матрицы по строкам 
    std::vector<int> colIndex; // Массив номеров столбцов 
    std::vector<int> rowPtr; // Массив индексов начала строк

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

// nz - n = четное
void InitCRSMatrix(CRSMatrix &matrix, int n, int nz)
{
    std::vector<double> initVec;
    const int varNum = 9;

    // 1. Сгенерировать вектор размерности ((nz - n)/2) без 0
    GenVecWithoutNull(initVec, (nz - n) / 2, varNum);

    // 2. Рандомно разместить его в верхнем треугольнике
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

    // 3. Заполнить Правильно главную диагональ
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
