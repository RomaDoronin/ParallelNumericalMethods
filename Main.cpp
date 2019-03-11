#include <iostream>
#include <vector>

using namespace std;

enum PNMStatus
{
	PNMStatusOk,
	PNMStatusFailed,
	PNMStatusInitError
};

// Блочное разложение Холецкого
// AMatrix - симетричная положительно определенная матрица
// LMatirx - нижнетреугольная с положительными элементами на диагонали
PNMStatus CholeskyBlockDecomposition(vector<int> AMatrix, size_t size, vector<int> &LMatirx)
{


	return PNMStatusOk;
}

void PrintMatrix(vector<int> vec, size_t size)
{
	size_t count = 0;

	for each (int val in vec)
	{
		cout << val << " ";
		if (count % size == 2)
		{
			cout << endl;
		}
		count++;
	}
}

bool CompareMatrix(vector<int> matrix1, vector<int> matrix2)
{
	for (int i = 0; i < matrix1.size(); i++)
	{
		if (matrix1[i] != matrix2[i])
			return false;
	}

	return true;
}

vector<int> GetTrMatrix(vector<int> matrix, size_t size)
{
	vector<int> trMatrix(matrix.size());

	for (int i = 0; i < matrix.size(); i++)
	{
		trMatrix[i / size + size * (i % size)] = matrix[i];
	}

	return trMatrix;
}

vector<int> GetMatrixRow(vector<int> matrix, size_t size, size_t rowNum)
{
	vector<int> row;

	for (int i = size * rowNum; i < size * (rowNum + 1); i++)
		row.push_back(matrix[i]);

	return row;
}

vector<int> GetMatrixColumn(vector<int> matrix, size_t size, size_t colNum)
{
	vector<int> col;

	for (int i = colNum; i < size * size; i += size)
		col.push_back(matrix[i]);

	return col;
}

vector<int> MultMatrix(vector<int> matrix1, vector<int> matrix2, size_t size)
{
	return vector<int>(2);
}

int main()
{
	vector<int> AMatrix = { 2,4,9, 4,8,0, 9,0,6 };
	size_t size = 3;
	vector<int> LMatirx;

	if (CholeskyBlockDecomposition(AMatrix, size, LMatirx) == PNMStatusOk)
	{
		cout << "LMatrix: " << endl;
		PrintMatrix(LMatirx, size);
	}
	else
	{
		cout << "Something goes wrong :(";
	}

	return 0;
}