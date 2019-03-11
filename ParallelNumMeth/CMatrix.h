#pragma once

#include <vector>
#include <iostream>

template <typename T>
class CMatrix
{
private:
    std::vector<T> m_matrix;
    size_t m_size;

#ifdef TEST_MODE
public:
#endif

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

        for each (T val in m_matrix)
        {
            std::cout << val << "	";
            if (count % m_size == 2)
            {
                std::cout << std::endl;
            }
            count++;
        }
    }

	CMatrix<T> GetTrMatrix() const
    {
        CMatrix<double> trMatrix(m_matrix.size());

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

	const bool operator == (const CMatrix& cMatrix)
	{
		for (int i = 0; i < m_size; i++)
		{
			if (m_matrix[i] != cMatrix.m_matrix[i])
			{
				return false;
			}
		}

		return true;
	}

};

