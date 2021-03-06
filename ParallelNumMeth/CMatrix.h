#pragma once

#ifndef _VECTOR_
#include <vector>
#endif
#ifndef _IOSTREAM_
#include <iostream>
#endif
#ifndef _FSTREAM_
#include <fstream>
#endif
#ifndef _STRING_
#include <string>
#endif
#ifndef _IOMANIP_
#include <iomanip>
#endif
#define _CMATRIX_


template <typename T>
class CMatrix
{
#ifdef TEST_MODE
public:
#else
private:
#endif
    std::vector<T> m_matrix;
    size_t m_size;

    std::vector<T> GetMatrixRow(size_t rowNum) const
    {
        std::vector<T> row;

        for (int i = m_size * rowNum; i < m_size * (rowNum + 1); i++)
            row.push_back(m_matrix[i]);

        return row;
    }

    std::vector<T> GetMatrixColumn(size_t colNum) const
    {
        std::vector<T> col;

        for (int i = colNum; i < m_size * m_size; i += m_size)
            col.push_back(m_matrix[i]);

        return col;
    }

    T ScalarMult(std::vector<T> vec1, std::vector<T> vec2)
    {
        T sum = 0;

        for (int i = 0; i < vec1.size(); i++)
        {
            sum += vec1[i] * vec2[i];
        }

        return sum;
    }

    void InitMatrix()
    {
        for (int i = 0; i < m_size * m_size; i++)
        {
            m_matrix.push_back(0);
        }
    }

    std::vector<T> GetSubMatrix(std::vector<T> matrix, size_t removeRow, size_t removeCol)
    {
        std::vector<T> submatrix;
        size_t size = sqrt(matrix.size());

        for (size_t i = 0; i < pow(size, 2); i++)
        {
            if ((i / size != removeRow) && (i % size != removeCol))
            {
                submatrix.push_back(matrix[i]);
            }
        }

        return submatrix;
    }

    T GetDeterminantRec(std::vector<T> matrix, size_t size)
    {
        if (size == 1)
            return matrix[0];

        T det = 0;

        for (size_t i = 0; i < size; i++)
        {
            det += matrix[i] * pow(-1, i) * GetDeterminantRec(GetSubMatrix(matrix, 0, i), size - 1);
        }

        return det;
    }

public:
    CMatrix(size_t size)
    {
        m_size = size;
        InitMatrix();
    }

    CMatrix(std::vector<T> matrix)
    {
        m_size = (size_t)sqrt(matrix.size());
        m_matrix = matrix;
    }

    CMatrix(const CMatrix<T> &matrix)
    {
        m_size = matrix.GetSize();

        for (int i = 0; i < pow(m_size, 2); i++)
        {
            m_matrix.push_back(matrix.m_matrix[i]);
        }
    }

    ~CMatrix()
    {
    }

    size_t GetSize() const
    {
        return m_size;
    }

    void SetMatrixCell(size_t index, T val)
    {
        m_matrix[index] = val;
    }

    void PrintMatrix() const
    {
        size_t count = 0;

        for (int i = 0; i < m_matrix.size(); i++)
        {
            std::cout << m_matrix[i] << "	";
            if (count % m_size == m_size - 1)
            {
                std::cout << std::endl;
            }
            count++;
        }
    }

    void WriteMatrixInFile(std::string filename, std::string subText) const
    {
        size_t count = 0;
        std::ofstream fout;
        fout.open(filename, std::ios_base::app);

        fout << "---------------------------------------- " << subText << std::endl;
        for (int i = 0; i < m_matrix.size(); i++)
        {
            fout << std::setprecision(3) << m_matrix[i] << "	";
            if (count % m_size == m_size - 1)
            {
                fout << std::endl;
            }
            count++;
        }

        fout.close();
    }

    T GetDeterminant()
    {
        return GetDeterminantRec(m_matrix, m_size);
    }

    bool IsPositivelyDefined()
    {
        std::vector<T> matrix = m_matrix;

        for (size_t size = m_size; size > 0; size--)
        {
            T det = GetDeterminantRec(matrix, size);
#ifdef TEST_MODE
            std::cout << "GetDeterminantRec(matrix, size) = " << det << std::endl;
#endif // TEST_MODE
            if (det <= 0)
                return false;

            matrix = GetSubMatrix(matrix, size - 1, size - 1);
        }

        return true;
    }

    CMatrix<T> GetTrMatrix() const
    {
        CMatrix<double> trMatrix(m_size);

        for (int i = 0; i < m_matrix.size(); i++)
        {
            trMatrix.SetMatrixCell(i / m_size + m_size * (i % m_size), m_matrix[i]);
        }

        return trMatrix;
    }

    T operator [] (int index)
    {
        return m_matrix[index];
    }

    CMatrix operator * (const CMatrix& cMatrix)
    {
        std::vector<T> matrix;

        for (int i = 0; i < m_matrix.size(); i++)
        {
            matrix.push_back(ScalarMult(GetMatrixRow(i / m_size), cMatrix.GetMatrixColumn(i % m_size)));
        }
        return CMatrix(matrix);
    }

    std::vector<T> operator * (const std::vector<T>& cVector)
    {
        std::vector<T> res;

        for (int i = 0; i < m_size; i++)
        {
            res.push_back(ScalarMult(cVector, GetMatrixRow(i)));
        }

        return res;
    }

    const bool operator == (const CMatrix& cMatrix)
    {
        for (int i = 0; i < pow(m_size, 2); i++)
        {
            if (m_matrix[i] != cMatrix.m_matrix[i])
            {
                return false;
            }
        }

        return true;
    }
};
