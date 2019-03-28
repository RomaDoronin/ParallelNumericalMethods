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

void GenVecWithoutNull(std::vector<double> &vec, int n, int var);

// nz - n = четное
void InitCRSMatrix(CRSMatrix &matrix, int n, int nz);

void PrintCRSMatrix(CRSMatrix &matrix);
