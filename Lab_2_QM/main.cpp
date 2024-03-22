#include <iostream>
#include <valarray>
#include <chrono>

using namespace std;
using real = double;
const int n = 1500;
const int m = 5;// номер  в группе
const int k = 1;// номер группы
//double A[n][n] ;// исходная матрица
//double X[n][1];//точное решение
real *A;// исходная матрица
real *X;//точное решение
real *B;

// метод для заполнения симметричной матрицы
void FillMatrixSymmetrically(real *&A, int n, int k) {
    A = new real[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
           // int random = rand() % 1001 - 1000;// Генерируем случайное число от -1000  до 0

            A[i * n + j] =  (double)rand() / RAND_MAX * -1000.0f;
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
            A[0] = sum + pow(10, 2 - k);
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
}

void A_to_LDLt(real *&A, int n) {
    real *D = new real[n];
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double d = A[j * n + i];
            D[j] = A[j * n + i];
            double a = A[j * n + i] / A[i * n + i];
            A[j * n + i] /= A[i * n + i];
            for (int k = i + 1; k <= j; ++k) {
                double aa = A[j * n + k] - A[j * n + i] * D[k];
                A[j * n + k] -= A[j * n + i] * D[k];
            }
        }
    }
}

// метод для заполнения точного решения
void FillX(real *&X, int n, int m) {
    X = new real[n];
    for (int i = 0; i < n; i++) {
        X[i] = m;
        m++;
    }
}

void MatrixMultiplication(real *&B, real *A,
                          int first_num_rows, int first_num_col,
                          real *X,
                          int second_num_rows, int second_num_col) {

    B = new real[first_num_rows /* *second_num_col*/];
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            double elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                double a = A[i * first_num_col + k];
                double x = X[k * second_num_col + j];
                elem += A[i * first_num_col + k] * X[k * second_num_col + j];
            }
            B[i/* *second_num_col+j*/] = elem;
        }
    }
    // return B;
}

real *FillY(real *&Y, real *A, real *B, int n) {
    Y = new real[n];
    double y = B[0];
    Y[0] = B[0];
    for (int i = 1; i < n; i++) {
        double b = B[i];
        real sum = B[i];
        for (int k = 0; k <= i - 1; k++) {
            if (i > k) {
                double yy = Y[k];
                double a = A[(i) * n + k];
                sum -= Y[k] * A[(i) * n + k];
            }
        }
        Y[i] = sum;
    }
    return Y;
}

real *FillZ(real *&Z, real *Y, real *A, int n) {
    Z = new real[n];
    for (int i = 0; i < n; i++) {
        double y = Y[i];
        double d = A[i * n + i];
        double z = Y[i] / A[i * n + i];
        Z[i] = Y[i] / A[i * n + i];
    }
    return Z;
}

real *FillX_find(real *&X, real *A, real *Z, int n) {
    X = new real[n];
    double zz = Z[n - 1];
    X[n - 1] = Z[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        double z = Z[i];
        real sum = z;
        for (int k = i + 1; k < n; k++) {
            if (k > i) {
                double a = A[k * n + i];
                double x = X[k];
                sum -= X[k] * A[k * n + i];
            }
        }
        //  double z=Z[i];
        X[i] = sum;
    }
}

// Находит номер строки ведущего по столбцу
int findRowOfLeadingByCol(real *A, int n, int colNum) {
    int row = colNum;
    int leading = A[colNum * n + colNum];
    for (int i = colNum + 1; i < n; i++) {
        real nextElem = A[i * n + colNum];
        if (abs(nextElem) > abs(leading)) {
            leading = nextElem;
            row = i;
        }
    }

    return row;
}

// Обратный ход, используется формула (9) на с.3
void goInReverseOrder(real *A_copy, real *B_copy, real *X, int n) {
    X[n - 1] = B_copy[n - 1] / A_copy[n * n - 1];
    for (int i = n - 2; i >= 0; i--) {
        // Находим сумму, от a_i_(i+1) до a_i_(n - 1), где i = 0..n-1
        real s = 0;
        for (int j = i + 1; j < n; j++) {
            s += A_copy[i * n + j] * X[j];
        }

        X[i] = 1 / A_copy[i * n + i] * (B_copy[i] - s);
    }
}

void Gauss_withLeading(real *A_copy, real *B_copy, int n, real *X) {

    // Прямой ход с выбором ведущего, используя формулы на с.14
    for (int k = 0; k < n - 1; k++) {
        int rowOfLeading = findRowOfLeadingByCol(A_copy, n, k);

        // Перемещаем строку с ведущим на нужную позицию в матрице
        for (int i = 0; i < n; i++) {
            swap(A_copy[rowOfLeading * n + i], A_copy[k * n + i]);
        }
        swap(B_copy[rowOfLeading], B_copy[k]);

        for (int i = k + 1; i < n; i++) {
            // Находим l_ik = a_ik / a_kk
            real l = A_copy[i * n + k] / A_copy[k * n + k];

            B_copy[i] -= l * B_copy[k];
            for (int j = k + 1; j < n; j++) {
                A_copy[i * n + j] -= l * A_copy[k * n + j];
            }
        }
    }

    // Обратный ход
    goInReverseOrder(A_copy, B_copy, X, n);
}

void Gauss_withoutLeading(real *A_copy, real *B_copy, int n, real *X) {

    // Прямой ход без выбора ведущего
    for (int k = 0; k < n - 1; k++) {
        for (int i = k + 1; i < n; i++) {
            // Находим l_ik = a_ik / a_kk
            real l = A_copy[i * n + k] / A_copy[k * n + k];

            B_copy[i] -= l * B_copy[k];
            for (int j = k + 1; j < n; j++) {
                A_copy[i * n + j] -= l * A_copy[k * n + j];
            }
        }
    }

    // Обратный ход
    goInReverseOrder(A_copy, B_copy, X, n);
}

// Евклидова норма
double EuclideanNorm(real *X, int n) {
    real sum = 0;
    for (int i = 0; i < n; i++) {
        sum += pow(X[i], 2);
    }
    return sqrt(sum);
}

// относительная погрешность
void RelativeError(double relativeError, real *X, real *X_calculated) {
    real *Distinction = new real[n];
    for (int i = 0; i < n; i++) {
        Distinction[i] = X[i] - X_calculated[i];
    }

    cout<<EuclideanNorm(X_calculated,n)<<"------"<<endl;
    relativeError = EuclideanNorm(Distinction, n) / EuclideanNorm(X, n);
}

int main() {
    FillMatrixSymmetrically(A, n, k);
    FillX(X, n, m);
// Находим b=AX
    MatrixMultiplication(B, A, n, n, X, n, 1);

    // Решаем по методу Гаусса с выбором ведущего_______________

    real *X_withLeading = new real[n];// полученный X в методе Гаусса с выбором ведущего

    real *A_copy = new real[n * n];
    for (int i = 0; i < n * n; i++) {
        A_copy[i] = A[i];
    }

    real *B_copy = new real[n];
    for (int i = 0; i < n; i++) {
        B_copy[i] = B[i];
    }
    //расставили границы для времени
    auto start_Gauss_withLeading = chrono::high_resolution_clock::now();
    Gauss_withLeading(A_copy, B_copy, n, X_withLeading);
    auto end_Gauss_withLeading = chrono::high_resolution_clock::now();
    int time_Gauss_withLeading = chrono::duration_cast<chrono::milliseconds>(
            end_Gauss_withLeading - start_Gauss_withLeading).count();
    // Конец решения по методу Гаусса с выбором ведущего______________


    // Решаем по методу Гаусса без выбора ведущего____________________

    real *X_withoutLeading = new real[n];// полученный X в методе Гаусса(без выбора ведущего)
    for (int i = 0; i < n * n; i++) {
        A_copy[i] = A[i];
    }

    for (int i = 0; i < n; i++) {
        B_copy[i] = B[i];
    }
    //расставили границы для времени
    auto start_Gauss_withoutLeading = chrono::high_resolution_clock::now();
    Gauss_withoutLeading(A_copy, B_copy, n, X_withoutLeading);
    auto end_Gauss_withoutLeading = chrono::high_resolution_clock::now();
    int time_Gauss_withoutLeading = chrono::duration_cast<chrono::milliseconds>(
            end_Gauss_withoutLeading - start_Gauss_withoutLeading).count();
    // Конец решения по методу Гаусса без выбора ведущего_____________

    // Решаем по методу LDLt_____________________________________

    real *Y = new real[n];
    real *Z = new real[n];
    real *Slay_X = new real[n];// полученный X в LDLt
    //расставили границы для времени
    auto start_LDLt = chrono::high_resolution_clock::now();
    A_to_LDLt(A, n);
    FillY(Y, A, B, n);
    FillZ(Z, Y, A, n);
    FillX_find(Slay_X, A, Z, n);
    auto end_LDLt = chrono::high_resolution_clock::now();
    int time_LDLt = chrono::duration_cast<chrono::milliseconds>(end_LDLt - start_LDLt).count();
    // Конец решения по методу LDLt_____________________________

    // Найдем относительные погрешности_____________________________
    double RelativeError_Gauss_withLeading;
    RelativeError(RelativeError_Gauss_withLeading, X, X_withLeading);

    double RelativeError_Gauss_withoutLeading;
    RelativeError(RelativeError_Gauss_withoutLeading, X, X_withoutLeading);

    double RelativeError_LDLt;
    RelativeError(RelativeError_LDLt, X, Slay_X);
    //_______________________________________________________________
    cout << "X: ";
    for (int i = 0; i < 5; i++) {
        cout << X[i] << ", ";
    }
    cout << endl << endl;
    // Вывод ответа в формате:
    //1. Первые 5 координат вектора приближённого решения x*.
    //2. Относительная погрешность
    //3. Время выполнения (можно приблизительно).


    cout << "The first 5 coordinates of the approximate solution vector x*:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << X_withLeading[i] << ", ";
    }
    cout << "    Gauss with leading" << endl;
    for (int i = 0; i < 5; i++) {
        cout << X_withLeading[i] << ", ";
    }
    cout << "    Gauss without leading" << endl;
    for (int i = 0; i < 5; i++) {
        cout << Slay_X[i] << ", ";
    }
    cout << "    LDLt" << endl;
    cout << endl;
    cout << "Relative error:" << endl;
    cout << RelativeError_Gauss_withLeading << "    Gauss with leading" << endl;
    cout << RelativeError_Gauss_withoutLeading << "    Gauss without leading" << endl;
    cout << RelativeError_LDLt << "    LDLt" << endl;
    cout << endl;
    cout << "Time:" << endl;
    cout << time_Gauss_withLeading << "    Gauss with leading" << endl;
    cout << time_Gauss_withoutLeading << "    Gauss without leading" << endl;
    cout << time_LDLt << "    LDLt" << endl;


}

