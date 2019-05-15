#include <gtest\gtest.h>

#define TEST_MODE

#include "..\ParallelNumMeth\CMatrix.h"

//---------------------------------------------------------------
// A22 - L21*L21T
double * WaveOperation(double * A22, double * L21, int n, int r)
{
    int size = n - r;
    double * res = new double[size * size];

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            int multVec = 0;
            for (int count = 0; count < r; count++)
            {
                //std::cout << "L21[i * r + count] : " << L21[i * r + count] << " | L21[j * r + count] : " << L21[j * r + count] << std::endl;
                multVec += L21[i * r + count] * L21[j * r + count];
            }
            //std::cout << "A22[i * size + j] : " << A22[i * size + j] << " | multVec : " << multVec << std::endl;
            res[i * size + j] = A22[i * size + j] - multVec;
        }
    }

    return res;
}

// Решается множество систем линейных уравнений: L21 * L11T = A21 -> L11 * L21T = A21T
void CulcL21(double * A21, double * L21, double * L11, int n, int r)
{
    int size = n - r;

    // Преведение всей системы к треугольному виду
    /*for (int mainDiagCount = 0; mainDiagCount < r; mainDiagCount++)
    {
        for (int rowCount = mainDiagCount + 1; rowCount < r; rowCount++)
        {
            int d = L11[rowCount * r + mainDiagCount] / L11[mainDiagCount * r + mainDiagCount];

            for (int colCount = mainDiagCount; colCount < r; colCount++)
            {
                L11[rowCount * r + colCount] -= d * L11[mainDiagCount * r + colCount];
            }

            for (int count = 0; count < size; count++)
            {
                //A21[rowCount * size + count] -= d * A21[mainDiagCount * r + count];
                // A21T
                A21[count * size + rowCount] -= d * A21[count * r + mainDiagCount];
            }
        }
    }*/

    // Обратный ход для каждой правой части
    for (int numRightPart = 0; numRightPart < size; numRightPart++)
    {
        for (int colCount = 0; colCount < r; colCount++)
        {
            int sum = 0;

            for (int rowCount = 0; rowCount < colCount; rowCount++)
            {
                // L11T
                sum += L11[colCount * r + rowCount] * L21[numRightPart * r + rowCount];
            }

            //std::cout << "A21[numRightPart * r + colCount] : " << A21[numRightPart * r + colCount] << " | sum : " << sum << " | L11[colCount * r + colCount] : " << L11[colCount * r + colCount] << std::endl;
            L21[numRightPart * r + colCount] = (A21[numRightPart * r + colCount] - sum) / L11[colCount * r + colCount];
        }
    }
}

void MergeLMatrix(double * L, double * L11, double * L21, double * L22, int n, int blockSize)
{
    int count = 0;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i < blockSize && j < blockSize)
            {
                L[count] = L11[i * blockSize + j];
            }
            else if (i < blockSize && j >= blockSize)
            {
                L[count] = 0;
            }
            else if (i >= blockSize && j < blockSize)
            {
                L[count] = L21[(i - blockSize) * blockSize + j];
            }
            else if (i >= blockSize && j >= blockSize)
            {
                L[count] = L22[(i - blockSize) * (n - blockSize) + (j - blockSize)];
            }

            count++;
        }
    }
}

bool PRKK(double * res1, double * res2, int n, double accuracy)
{
    for (int i = 0; i < n; i++)
    {
        if (abs(res1[i] - res2[i]) > accuracy)
            return false;
    }

    return true;
}
//---------------------------------------------------------------

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

//------------------- WaveOperation
TEST(Functions, WaveOperation)
{
    int n, r;
    double * A22;
    double * L21;
    double * res;
    double * refRes;

    // --- Test 1
    n = 4;
    r = 2;
    A22 = new double[(n - r) * (n - r)] { 2,3,
                                          4,5 };
    L21 = new double[(n - r) * r] { 6,0,
                                    3,4 };

    res = WaveOperation(A22, L21, n, r);
    refRes = new double[n - r] { -34,-15, -14,-20 };

    EXPECT_EQ(PRKK(res, refRes, (n - r) * (n - r), 0.0001), true);

    // --- Test 2
    n = 9;
    r = 3;
    A22 = new double[(n - r) * (n - r)]{ 100,100,100,100,100,100,
                                         100,100,100,100,100,100,
                                         100,100,100,100,100,100,
                                         100,100,100,100,100,100,
                                         100,100,100,100,100,100,
                                         100,100,100,100,100,100, };
    L21 = new double[(n - r) * r]{ 6,0,3,
                                   4,7,5,
                                   3,6,5,
                                   8,4,8,
                                   8,4,7,
                                   5,6,2 };

    res = WaveOperation(A22, L21, n, r);
    refRes = new double[n - r]{ 55,61,67,28,31,64,
                                61,10,21,0,5,28,
                                67,21,30,12,17,39,
                                28,0,12,-44,-36,20,
                                31,5,17,-36,-29,22,
                                64,28,39,20,22,35 };

    EXPECT_EQ(PRKK(res, refRes, (n - r) * (n - r), 0.0001), true);
}

//------------------- CulcL21
TEST(Functions, CulcL21)
{
    int n, r;
    double * L21;
    double * L11;
    double * A21;
    double * refL21;

    // --- Test 1
    n = 5;
    r = 2;

    L21 = new double[(n - r) * r];
    L11 = new double[r * r] { 1,0,
                              6,2 };
    A21 = new double[(n - r) * r] { 5,7,
                                    1,2,
                                    5,6};

    CulcL21(A21, L21, L11, n, r);

    refL21 = new double[(n - r) * r] { 5,-11.5,
                                       1,   -2,
                                       5,  -12 };

    EXPECT_EQ(PRKK(L21, refL21, (n - r) * r, 0.0001), true);

    // --- Test 2
    n = 9;
    r = 3;

    L21 = new double[(n - r) * r];
    L11 = new double[r * r]{ 1,0,0,
                             6,2,0,
                             8,4,7 };
    A21 = new double[(n - r) * r]{ 5,7,8,
                                   1,2,7,
                                   5,6,5,
                                   8,5,7,
                                   3,6,2,
                                   1,4,5 };

    CulcL21(A21, L21, L11, n, r);

    refL21 = new double[(n - r) * r]{ 5,-11.5,    2,
                                      1,   -2,    1,
                                      5,  -12,13./7.,
                                      8,-21.5,29./7.,
                                      3,   -6, 2./7.,
                                      1, -1, 1./7. };

    EXPECT_EQ(PRKK(L21, refL21, (n - r) * r, 0.0001), true);
}

//------------------- MergeLMatrix
TEST(Functions, MergeLMatrix)
{
    int n = 6;
    int r = 2;

    double * L = new double[n * n];
    double * L11 = new double[r * r] { 8,0,
                                       6,7 };

    double * L21 = new double[(n - r) * r] { 6,3,
                                             6,8,
                                             9,8,
                                             3,5, };

    double * L22 = new double[(n - r) * (n - r)]{ 6,0,0,0,
                                                  6,8,0,0,
                                                  9,8,5,0,
                                                  3,5,1,4 };

    MergeLMatrix(L, L11, L21, L22, n, r);

    double * refL = new double[n * n] { 8,0,0,0,0,0,
                                        6,7,0,0,0,0,
                                        6,3,6,0,0,0,
                                        6,8,6,8,0,0,
                                        9,8,9,8,5,0,
                                        3,5,3,5,1,4 };

    EXPECT_EQ(PRKK(L, refL, n * n, 0.0001), true);
}

int main(int argc, char* argv[])
{
    testing::InitGoogleTest(&argc, argv);
    int res = RUN_ALL_TESTS();
    int a;
    //std::cin >> a;
    return res;
}
