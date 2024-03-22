//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "cert-msc50-cpp"

#include <iostream>
#include <valarray>
#include <chrono>

using namespace std;
using real = double;
const int n_ = 1500;
const int m_ = 5;// номер  в группе
const int k_ = 1;// номер группы
real *A_;// исходная матрица
real *Y_;//точное решение
real *F_;

void FillA(real *&A, int n, int m, int k) {
    A = new real[(n + 1) * (n + 1)];
    for(int i=0;i<(n+1)*(n+1);i++){
        A[i]=0;
    }
    A[0] = m;
    //заполнили диагональ
    for (int i = 1; i < n + 1; i++) {
        A[i * (n + 1) + i] = m + k + i - 1;
    }
    //заполнили верхнюю диагональ
    for (int i = 0; i < n; i++) {
        A[i * (n + 1) + i + 1] = m + i - 1;
    }
    //заполнили нижнюю диагональ
    for (int i = 0; i < n; i++) {
        A[(i + 1) * (n + 1) + i] = -k;
    }
}

void FillY(real *&Y, int n) {
    Y = new real[n + 1];
    for (int i = 0; i < n + 1; i++) {
        Y[i] = i + 1;
    }
}

void MatrixMultiplication(real *&F, const real *A,
                          int first_num_rows, int first_num_col,
                          const real *X,
                          int second_num_rows, int second_num_col) {

    F = new real[first_num_rows];
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                elem += A[i * first_num_col + k] * X[k * second_num_col + j];
            }
            F[i] = elem;
        }
    }
}


// Находит номер строки ведущего по столбцу
int findRowOfLeadingByCol(const real *A, int n, int colNum) {
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
void goInReverseOrder(const real *A_copy, const real *B_copy, real *X, int n) {
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
double RelativeError(double relativeError, real *X, const real *X_calculated, int n) {
    real *Distinction = new real[n];
    for (int i = 0; i < n; i++) {
        Distinction[i] = X[i] - X_calculated[i];
    }
    relativeError = EuclideanNorm(Distinction, n) / EuclideanNorm(X, n);
    return relativeError;
}




real *Find_Slay_Y(real *&Y, const real *A, const real *F, int n) {
    real *Alf = new real[n + 1];
    real *B = new real[n + 1];
    Y = new real[n + 1];
    //прямая прогонка
    Alf[0] = -A[1] / A[0];
    B[0] = F[0] / A[0];
    real denominator;
    for (int i = 1; i < n; i++) {
        denominator = A[i * (n + 1) + i] - (-A[i * (n + 1) + (i - 1)]) * Alf[i - 1];
        Alf[i] = -A[i * (n + 1) + i + 1] / denominator;
        B[i] = (F[i] + (-A[i * (n + 1) + (i - 1)]) * B[i - 1]) / denominator;
    }
    B[n] = (F[n] + (-A[(n + 1) * (n + 1) - 2]) * B[n - 1]) /
           (A[(n + 1) * (n + 1) - 1] - (-A[(n + 1) * (n + 1) - 2]) * Alf[n - 1]);

    //обратная прогонка
    Y[n] = B[n];
    for (int i = n - 1; i >= 0; i--) {
        Y[i] = Alf[i] * Y[i + 1] + B[i];
    }
    return Y;
}

int main() {
    FillA(A_, n_, m_, k_);
    FillY(Y_, n_);
    // Находим f=Ay
    MatrixMultiplication(F_, A_, n_ + 1, n_ + 1, Y_, n_ + 1, 1);

    // Решаем по методу Гаусса с выбором ведущего_______________

    real *Y_withLeading = new real[n_+1];// полученный Y в методе Гаусса с выбором ведущего

    real *A_copy = new real[(n_+1) * (n_+1)];
    for (int i = 0; i < (n_+1) * (n_+1); i++) {
        A_copy[i] = A_[i];
    }

    real *F_copy = new real[n_+1];
    for (int i = 0; i < n_+1; i++) {
        F_copy[i] = F_[i];
    }
    // расставили границы для времени
    auto start_Gauss_withLeading = chrono::high_resolution_clock::now();
    Gauss_withLeading(A_copy, F_copy, n_+1, Y_withLeading);
    auto end_Gauss_withLeading = chrono::high_resolution_clock::now();
    int time_Gauss_withLeading = chrono::duration_cast<chrono::milliseconds>(
            end_Gauss_withLeading - start_Gauss_withLeading).count();
    //  Конец решения по методу Гаусса с выбором ведущего______________


    //  Решаем по методу Гаусса без выбора ведущего____________________

    real *Y_withoutLeading = new real[n_+1];// полученный Y в методе Гаусса(без выбора ведущего)
    for (int i = 0; i < (n_+1) * (n_+1); i++) {
        A_copy[i] = A_[i];
    }

    for (int i = 0; i < n_+1; i++) {
        F_copy[i] = F_[i];
    }
    // расставили границы для времени
    auto start_Gauss_withoutLeading = chrono::high_resolution_clock::now();
    Gauss_withoutLeading(A_copy, F_copy, n_+1, Y_withoutLeading);
    auto end_Gauss_withoutLeading = chrono::high_resolution_clock::now();
    int time_Gauss_withoutLeading = chrono::duration_cast<chrono::milliseconds>(
            end_Gauss_withoutLeading - start_Gauss_withoutLeading).count();
    //  Конец решения по методу Гаусса без выбора ведущего_____________


    //   Решаем по методу прогонки_____________________________________

    real *Slay_Y = new real[n_ + 1];// полученный Y в методе прогонки
    //расставили границы для времени
    auto start_Sweep_method = chrono::high_resolution_clock::now();
    Find_Slay_Y(Slay_Y, A_, F_, n_);
    auto end_Sweep_method = chrono::high_resolution_clock::now();
    int time_Sweep_method = chrono::duration_cast<chrono::milliseconds>(end_Sweep_method - start_Sweep_method).count();
    // Конец решения по методу прогонки_____________________________
    // Найдем относительные погрешности_____________________________
    double RelativeError_Gauss_withLeading = 0;
    RelativeError_Gauss_withLeading=RelativeError(RelativeError_Gauss_withLeading, Y_, Y_withLeading,n_+1);

    double RelativeError_Gauss_withoutLeading=0;
    RelativeError_Gauss_withoutLeading=RelativeError(RelativeError_Gauss_withoutLeading, Y_, Y_withoutLeading,n_+1);

    double RelativeError_Sweep_method=0;
    RelativeError_Sweep_method=RelativeError(RelativeError_Sweep_method,Y_, Slay_Y,n_+1);
    //_______________________________________________________________
    cout << "Y: ";
    for (int i = 0; i < 5; i++) {
        cout << Y_[i] << ", ";
    }
    cout << endl << endl;


    // Вывод ответа в формате:
    //1. Первые 5 координат вектора приближённого решения x*.
    //2. Относительная погрешность
    //3. Время выполнения (можно приблизительно).


    cout << "The first 5 coordinates of the approximate solution vector x*:" << endl;
    for (int i = 0; i < 5; i++) {
        cout << Y_withLeading[i] << ", ";
    }
    cout << "    Gauss with leading" << endl;
    for (int i = 0; i < 5; i++) {
        cout << Y_withLeading[i] << ", ";
    }
    cout << "    Gauss without leading" << endl;
    for (int i = 0; i < 5; i++) {
        cout << Slay_Y[i] << ", ";
    }
    cout << "    Sweep method" << endl;


    cout << endl;
    cout << "Relative error:" << endl;
    cout << RelativeError_Gauss_withLeading << "    Gauss with leading" << endl;
    cout << RelativeError_Gauss_withoutLeading << "    Gauss without leading" << endl;
    cout << RelativeError_Sweep_method << "    Sweep method" << endl;
    cout << endl;
    cout << "Time:" << endl;
    cout << time_Gauss_withLeading << "    Gauss with leading" << endl;
    cout << time_Gauss_withoutLeading << "    Gauss without leading" << endl;
    cout << time_Sweep_method << "    Sweep method" << endl;
}