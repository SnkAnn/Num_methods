#include <chrono>
#include <iostream>
#include "stdlib.h"
#include <math.h>

using namespace std;

//using real = float;
using real = double;
real *A; // Матрица А
const int n = 1500; // Порядок матрицы А
const int m = 5;// Номер в списке группы



real *B; // Правая часть уравнения AX = B

real *X_precise; // Точное решение AX = B

real* multiplyMatrices(real *first, int firstRows, int firstCols,
                       real *second, int secondRows, int secondCols)
{
    real * result = new real[firstRows * secondCols];
    for (int i = 0; i < firstRows; i++)
    {
        for (int k = 0; k < secondCols; k++)
        {
            real elem = 0;
            for (int j = 0; j < firstCols; j++)
            {
                elem += first[i * firstCols + j] * second[j * secondCols + k];
            }
            result[i * secondCols + k] = elem;
        }
    }

    return result;
}

real* subtractMatrix(real *from, real *toSubtract, int n, int m)
{
    real* result = new real[n * m];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            int idx = i * m + j;
            result[idx] = from[idx] - toSubtract[idx];
        }
    }

    return result;
}

// Находит невязку
real* findResidual(real *A, real *X, real *B, int n)
{

    real* B_approximate = multiplyMatrices(A, n, n,
                                           X, n, 1);

    real* residual = subtractMatrix(B, B_approximate, n, 1);

    delete[] B_approximate;

    return residual;
}

real findMaximumNorm(real *matrix, int n, int m)
{
    real norm = 0;
    for (int j = 0; j < n; j++)
    {
        real sum = 0;
        for (int i = 0; i < m; i++)
        {
            sum += abs(matrix[j * m + i]);
        }

        if (sum > norm)
        {
            norm = sum;
        }
    }

    return norm;
}

// Находит приблизительную погрешность
real findRelativeError(real* X_approximate, real* X_precise, int n) {
    real *diff = subtractMatrix(X_precise, X_approximate, n, 1);

    // Относительная погрешность
    real relError = findMaximumNorm(diff, n, 1) / findMaximumNorm(X_precise, n, 1);

    delete[] diff;

    return relError;
}

// Заполняет исходные матрицы (Lab 1)
void initValuesForLab_1(real *&A, int n, real *&B, real *&X_precise, int m) {
        A = new real[n * n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                int random = rand() % 1001 - 1000;// Генерируем случайное число от -1000  до 0

                A[i * n + j] = random;
            }
        }
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                A[i * n + j] = A[j * n + i];
            }
        }
        for (int i = 0; i < n; i++) {
            if (i == 0) {
                real sum = 0;
                int j = 1;
                while (j != n) {
                    double a = A[i * n + j];
                    sum -= A[i * n + j];
                    j++;
                }
                A[0] = sum + pow(10, 2 - 1);
            } else {
                real sum = 0;
                int j = 0;
                while (j != n) {
                    if (j != i) {
                        double a = A[i * n + j];
                        sum -= A[i * n + j];
                    }
                    j++;
                }
                A[i * n + i] = sum;
            }
        }

    X_precise = new real[n];
    for (int i = m; i < n + m; i++)
    {
        X_precise[i - m] = i; // Заполняем точное решение
    }

    B = multiplyMatrices(A, n, n,
                         X_precise, n, 1);
}

// Заполняет исходные матрицы (Лаб. №2)
void initValuesForLab_2(real *&A, int n, real *&B, real *&X_precise, int groupNumber, int m)
{
    A = new real[n * n];
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)  // Заполняем лишь нижний треугольник
                break;
            int random = rand() % 1001 - 1000; // Генерируем случайное число от 0  до -1000
            A[i * n + j] = random;
        }
    }

    // Для удобства заполняем верхний треугольник
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            A[i * n + j] = A[j * n + i];
        }
    }

    // Инициализируем диагональные
    real a00 = 0;
    for (int i = 1; i < n; i++) {
        a00 -= A[i];
    }
    a00 += pow(10, 2 - groupNumber);
    A[0] = a00;

    for (int i = 1; i < n; i++) {
        real a_ii = 0;
        for (int j = 0; j < n; j++) {
            if (j == i)
                continue;
            a_ii -= A[i * n + j];
        }
        A[i * n + i] = a_ii;
    }

    X_precise = new real[n];

    // Заполняем точное решение
    for (int i = m; i < n + m; i++)
    {
        X_precise[i - m] = i;
    }

    B = multiplyMatrices(A, n, n,
                         X_precise, n, 1);
}

// Находит номер строки ведущего по столбцу
int findRowOfLeadingByCol(real *A, int n, int colNum)
{
    int row = colNum;
    int leading = A[colNum * n + colNum];
    for (int i = colNum + 1; i < n; i++)
    {
        real nextElem = A[i * n + colNum];
        if (abs(nextElem) > abs(leading))
        {
            leading = nextElem;
            row = i;
        }
    }

    return row;
}

// Обратный ход, используется формула (9) на с.3
void goInReverseOrder(real *A_copy, real *B_copy, real *X, int n)
{
    X[n - 1] = B_copy[n - 1] / A_copy[n * n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        // Находим сумму, от a_i_(i+1) до a_i_(n - 1), где i = 0..n-1
        real s = 0;
        for (int j = i + 1; j < n; j++)
        {
            s += A_copy[i * n + j] * X[j];
        }

        X[i] = 1 / A_copy[i * n + i] * (B_copy[i] - s);
    }
}

// Используются формулы (6)-(9) на с.3
// Возвращает время, затраченное на вычисления (millisec)
int solveUsingGaussNoLeading(real *A, real *B, int n, real *X)
{
    real *A_copy = new real[n * n];
    for (int i = 0; i < n * n; i++)
    {
        A_copy[i] = A[i];
    }

    real *B_copy = new real[n];
    for (int i = 0; i < n; i++)
    {
        B_copy[i] = B[i];
    }

    auto start = chrono::high_resolution_clock::now();

    // Прямой ход без выбора ведущего
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            // Находим l_ik = a_ik / a_kk
            real l = A_copy[i * n + k] / A_copy[k * n + k];

            B_copy[i] -= l * B_copy[k];
            for (int j = k + 1; j < n; j++)
            {
                A_copy[i * n + j] -= l * A_copy[k * n + j];
            }
        }
    }

    // Обратный ход
    goInReverseOrder(A_copy, B_copy, X, n);
    auto end = chrono::high_resolution_clock::now();

    delete[] A_copy;
    delete[] B_copy;

    return chrono::duration_cast<chrono::milliseconds>(end - start).count();
}

// Возвращает время, затраченное на вычисления (millisec)
int solveUsingGaussWithLeading(real *A, real *B, int n, real *X)
{
    real *A_copy = new real[n * n];
    for (int i = 0; i < n * n; i++)
    {
        A_copy[i] = A[i];
    }

    real *B_copy = new real[n];
    for (int i = 0; i < n; i++)
    {
        B_copy[i] = B[i];
    }

    auto start = chrono::high_resolution_clock::now();

    // Прямой ход с выбором ведущего, используя формулы на с.14
    for (int k = 0; k < n - 1; k++)
    {
        int rowOfLeading = findRowOfLeadingByCol(A_copy, n, k);

        // Перемещаем строку с ведущим на нужную позицию в матрице
        for (int i = 0; i < n; i++)
        {
            swap(A_copy[rowOfLeading * n + i], A_copy[k * n + i]);
        }
        swap(B_copy[rowOfLeading], B_copy[k]);

        for (int i = k + 1; i < n; i++)
        {
            // Находим l_ik = a_ik / a_kk
            real l = A_copy[i * n + k] / A_copy[k * n + k];

            B_copy[i] -= l * B_copy[k];
            for (int j = k + 1; j < n; j++)
            {
                A_copy[i * n + j] -= l * A_copy[k * n + j];
            }
        }
    }

    // Обратный ход
    goInReverseOrder(A_copy, B_copy, X, n);
    auto end = chrono::high_resolution_clock::now();

    delete[] A_copy;
    delete[] B_copy;

    return chrono::duration_cast<chrono::milliseconds>(end - start).count();
}

int main()
{
    initValuesForLab_1(A, n, B, X_precise, m);
    // initValuesForLab_2(A, n, B, X_precise, 1, m);

    // Решения, получившиеся после метода Гаусса
    real *X_withLeading = new real[n];
    real *X_noLeading = new real[n];

    int timeElapsedWithLeading = solveUsingGaussWithLeading(A, B, n, X_withLeading);
    int timeElapsedNoLeading = solveUsingGaussNoLeading(A, B, n, X_noLeading);

    // Невязки
    real *residualWithLeading = findResidual(A, X_withLeading, B, n);
    real *residualNoLeading = findResidual(A, X_noLeading, B, n);

    // Относительная погрешность с ведущим
    real relErrorWithLeading = findRelativeError(X_withLeading, X_precise, n);
    // Относительная погрешность без ведущего
    real relErrorNoLeading = findRelativeError(X_noLeading, X_precise, n);


    cout << "With the host: \n\tx*: (";
    for (int i = 0; i < 5; i++)
    {
        cout << X_withLeading[i] << ", ";
    }
    cout << "...)\n";
    cout << "\tMaximum norm of discrepancy: " << findMaximumNorm(residualWithLeading, n, 1) << "\n";
    cout << "\tRelative error: " << relErrorWithLeading << "\n";
    cout << "\tTime spent on calculations: " << timeElapsedWithLeading << '\n';

    //////////////////////////////////////
    cout << "Without a presenter: \n\tx*: (";
    for (int i = 0; i < 5; i++)
    {
        cout << X_noLeading[i] << ", ";
    }
    cout << "...)\n";
    cout << "\tMaximum norm of discrepancy: "  <<  findMaximumNorm(residualNoLeading, n, 1) << "\n";
    cout << "\tRelative error: "  << relErrorNoLeading << "\n";
    cout << "\t Time spent on calculations: " << timeElapsedNoLeading << '\n';
}