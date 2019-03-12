#include <gtest\gtest.h>

#define TEST_MODE

#include "..\ParallelNumMeth\CMatrix.h"

//-------------------- GetMatrixRow
TEST(Functions, GetMatrixRow)
{
	CMatrix<double> cMatrix2({ 5,4, 7,5 });
	EXPECT_EQ(cMatrix2.GetMatrixRow(0)[0], 5);
	EXPECT_EQ(cMatrix2.GetMatrixRow(0)[1], 4);
	EXPECT_EQ(cMatrix2.GetMatrixRow(1)[0], 7);
	EXPECT_EQ(cMatrix2.GetMatrixRow(1)[1], 5);

	CMatrix<double> cMatrix3({ 5,4,7, 5,3,9, 7,2,1 });
	EXPECT_EQ(cMatrix3.GetMatrixRow(0)[0], 5);
	EXPECT_EQ(cMatrix3.GetMatrixRow(0)[1], 4);
	EXPECT_EQ(cMatrix3.GetMatrixRow(0)[2], 7);
	EXPECT_EQ(cMatrix3.GetMatrixRow(1)[0], 5);
	EXPECT_EQ(cMatrix3.GetMatrixRow(1)[1], 3);
	EXPECT_EQ(cMatrix3.GetMatrixRow(1)[2], 9);
	EXPECT_EQ(cMatrix3.GetMatrixRow(2)[0], 7);
	EXPECT_EQ(cMatrix3.GetMatrixRow(2)[1], 2);
	EXPECT_EQ(cMatrix3.GetMatrixRow(2)[2], 1);
}

//-------------------- GetMatrixColumn
TEST(Functions, GetMatrixColumn)
{
	CMatrix<double> cMatrix2({ 5,4, 7,5 });
	EXPECT_EQ(cMatrix2.GetMatrixColumn(0)[0], 5);
	EXPECT_EQ(cMatrix2.GetMatrixColumn(0)[1], 7);
	EXPECT_EQ(cMatrix2.GetMatrixColumn(1)[0], 4);
	EXPECT_EQ(cMatrix2.GetMatrixColumn(1)[1], 5);

	CMatrix<double> cMatrix3({ 5,4,7, 5,3,9, 7,2,1 });
	EXPECT_EQ(cMatrix3.GetMatrixColumn(0)[0], 5);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(0)[1], 5);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(0)[2], 7);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(1)[0], 4);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(1)[1], 3);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(1)[2], 2);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(2)[0], 7);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(2)[1], 9);
	EXPECT_EQ(cMatrix3.GetMatrixColumn(2)[2], 1);
}

//-------------------- ScalarMult
TEST(Functions, ScalarMult)
{
	CMatrix<double> cMatrix(3);
	
	EXPECT_EQ(cMatrix.ScalarMult({ 4,3 }, { 2,6 }), 26);
	EXPECT_EQ(cMatrix.ScalarMult({ 4,2,3 }, { 2,0,6 }), 26);
	EXPECT_EQ(cMatrix.ScalarMult({ 4,2,3,-2 }, { 2,0,6,5 }), 16);
	EXPECT_EQ(cMatrix.ScalarMult({ -4,2,3,5,-3 }, { 2,0,6,-3,1 }), -8);
}

//-------------------- GetSize
TEST(Functions, GetSize)
{
	CMatrix<double> cMatrix1(5);
	EXPECT_EQ(cMatrix1.GetSize(), 5);
	CMatrix<double> cMatrix2({ 5,5, 63,-25 });
	EXPECT_EQ(cMatrix2.GetSize(), 2);
	CMatrix<double> cMatrix3({ 5,3,4, 7,5,-2, 5,-2,10 });
	EXPECT_EQ(cMatrix3.GetSize(), 3);
	CMatrix<double> cMatrix4({ 5,3,-6,4, 2,3,-6,10, 3,0,0,0, -25,2,-2,0 });
	EXPECT_EQ(cMatrix4.GetSize(), 4);
}

//-------------------- SetMatrixCell
TEST(Functions, SetMatrixCell)
{
	CMatrix<double> cMatrix2({ 5,4, 7,5 });
	cMatrix2.SetMatrixCell(0, 11);
	EXPECT_EQ(cMatrix2[0], 11);
	cMatrix2.SetMatrixCell(1, 11);
	EXPECT_EQ(cMatrix2[1], 11);
	cMatrix2.SetMatrixCell(2, 11);
	EXPECT_EQ(cMatrix2[2], 11);
	cMatrix2.SetMatrixCell(3, 11);
	EXPECT_EQ(cMatrix2[3], 11);

	CMatrix<double> cMatrix3({ 5,4,7, 5,3,9, 7,2,1 });
	cMatrix3.SetMatrixCell(0, 11);
	EXPECT_EQ(cMatrix3[0], 11);
	cMatrix3.SetMatrixCell(1, 11);
	EXPECT_EQ(cMatrix3[1], 11);
	cMatrix3.SetMatrixCell(2, 11);
	EXPECT_EQ(cMatrix3[2], 11);
	cMatrix3.SetMatrixCell(3, 11);
	EXPECT_EQ(cMatrix3[3], 11);
	cMatrix3.SetMatrixCell(4, 11);
	EXPECT_EQ(cMatrix3[4], 11);
	cMatrix3.SetMatrixCell(5, 11);
	EXPECT_EQ(cMatrix3[5], 11);
	cMatrix3.SetMatrixCell(6, 11);
	EXPECT_EQ(cMatrix3[6], 11);
	cMatrix3.SetMatrixCell(7, 11);
	EXPECT_EQ(cMatrix3[7], 11);
	cMatrix3.SetMatrixCell(8, 11);
	EXPECT_EQ(cMatrix3[8], 11);
}

//-------------------- GetTrMatrix
TEST(Functions, GetTrMatrix)
{
	CMatrix<double> cMatrix1({ 5,4, 7,5 });
	CMatrix<double> cMatrix2({ 5,7, 4,5 });
	EXPECT_EQ(cMatrix1.GetTrMatrix() == cMatrix2, true);

	CMatrix<double> cMatrix3({ 5,4,7, 5,3,9, 7,2,1 });
	CMatrix<double> cMatrix4({ 5,5,7, 4,3,2, 7,9,1 });
	EXPECT_EQ(cMatrix3.GetTrMatrix() == cMatrix4, true);
}

//-------------------- []
TEST(Functions, RightBrackets)
{
	CMatrix<double> cMatrix2({ 5,4, 7,5 });
	EXPECT_EQ(cMatrix2[0], 5);
	EXPECT_EQ(cMatrix2[1], 4);
	EXPECT_EQ(cMatrix2[2], 7);
	EXPECT_EQ(cMatrix2[3], 5);

	CMatrix<double> cMatrix3({ 5,4,7, 5,3,9, 7,2,1 });
	EXPECT_EQ(cMatrix3[0], 5);
	EXPECT_EQ(cMatrix3[1], 4);
	EXPECT_EQ(cMatrix3[2], 7);
	EXPECT_EQ(cMatrix3[3], 5);
	EXPECT_EQ(cMatrix3[4], 3);
	EXPECT_EQ(cMatrix3[5], 9);
	EXPECT_EQ(cMatrix3[6], 7);
	EXPECT_EQ(cMatrix3[7], 2);
	EXPECT_EQ(cMatrix3[8], 1);
}

//-------------------- *
TEST(Functions, Mult)
{
	CMatrix<double> cMatrix1({ 5,4, 7,5 });
	CMatrix<double> cMatrix2({ 1,0, 0,1 });
	CMatrix<double> cMatrix12({ 5,4, 7,5 });
	EXPECT_EQ((cMatrix1 * cMatrix2) == cMatrix12, true);

	CMatrix<double> cMatrix3({ 5,4, 7,5 });
	CMatrix<double> cMatrix4({ 0,0, 0,0 });
	CMatrix<double> cMatrix34({ 0,0, 0,0 });
	EXPECT_EQ((cMatrix3 * cMatrix4) == cMatrix34, true);

	CMatrix<double> cMatrix5({ 5,2,6, 8,4,7, 5,2,6 });
	CMatrix<double> cMatrix6({ 0,3,2, 1,4,5, 8,7,9 });
	CMatrix<double> cMatrix56({ 50,65,74, 60,89,99, 50,65,74 });
	EXPECT_EQ((cMatrix5 * cMatrix6) == cMatrix56, true);
}

//------------------- GetSubMatrix
TEST(Functions, GetSubMatrix)
{
	CMatrix<double> cMatrix1({ 5,2,6, 8,4,7, 5,2,6 });
	CMatrix<double> cMatrix2({ 5,2, 8,4 });
	CMatrix<double> cMatrix3(cMatrix1.GetSubMatrix(cMatrix1.m_matrix, 2, 2));
	EXPECT_EQ(cMatrix3 == cMatrix2, true);

	CMatrix<double> cMatrix4({ 5,2,6, 8,4,7, 5,2,6 });
	CMatrix<double> cMatrix5({ 5,6, 5,6 });
	CMatrix<double> cMatrix6(cMatrix4.GetSubMatrix(cMatrix4.m_matrix, 1, 1));
	EXPECT_EQ(cMatrix5 == cMatrix6, true);
}

//------------------- GetDeterminantRec
TEST(Functions, GetDeterminantRec)
{
	CMatrix<double> cMatrix1({ 2,4,9, 4,8,1, 9,1,6 });
	EXPECT_EQ(cMatrix1.GetDeterminant(), -578);

	std::vector<double> vec1 = { 2,4, 4,8 };
	EXPECT_EQ(cMatrix1.GetDeterminantRec(vec1, 2), 0);

	std::vector<double> vec2;
	vec2.push_back(2);
	EXPECT_EQ(cMatrix1.GetDeterminantRec(vec2, 1), 2);
}

//------------------- IsPositivelyDefined
TEST(Functions, IsPositivelyDefined)
{
	CMatrix<double> cMatrix1({ 2,4,9, 4,8,1, 9,1,6 });
	EXPECT_EQ(cMatrix1.IsPositivelyDefined(), false);

	CMatrix<double> cMatrix2({ 2,-1,0, -1,2,-1, 0,-1,2 });
	EXPECT_EQ(cMatrix2.IsPositivelyDefined(), true);

	CMatrix<double> cMatrix3({ 10,-3,2, -3,3,-2, 2,-2,7 });
	EXPECT_EQ(cMatrix3.IsPositivelyDefined(), true);
}

int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	int res = RUN_ALL_TESTS();
	int a;
	std::cin >> a;
	return res;
}