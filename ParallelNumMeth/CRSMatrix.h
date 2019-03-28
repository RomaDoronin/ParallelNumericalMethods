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

void GenVecWithoutNull(std::vector<double> &vec, int n, int var);

// nz - n = ������
void InitCRSMatrix(CRSMatrix &matrix, int n, int nz);

void PrintCRSMatrix(CRSMatrix &matrix);
