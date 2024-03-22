//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "cert-msc50-cpp"

#include <iostream>
#include <valarray>
#include <chrono>

using namespace std;
using real = double;
const int n_ = 12;
const int m_ = 5;// номер  в группе
real *A_;// исходная матрица
real *X_;//точное решение
real *F_;

void FillA(real *&A, int n) {
    A = new real[n * n];
    for (int i = 0; i < n; i++) {
        real sum = 0;
        for (int j = 0; j < n; j++) {
            if (i != j) {
                int random = rand() % 6 - 5;
                A[i * n + j] = random;
                sum += random;
            }
        }
        if (i == 0) {
            A[i * n + i] = -sum + 1;
        } else {
            A[i * n + i] = -sum;
        }
    }
}

void FillX(real *&X, int n, real m) {
    X = new real[n];
    for (int i = 0; i < n; i++) {
        X[i] = m + i;
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

real *Jacobi_met(real *A, real *F, real *&X1, real *&X2, int n, int k) {
    if (k % 2 == 1) {//нечетная итерация x1 прошлое x2
        real fun = 0;
        real sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i * n + j] * X1[j];
                }
            }
            fun = F[i] - sum;
            X2[i] = fun / A[i * n + i];
            sum = 0;
        }
        return X2;
    } else {//четная итерация x2 прошлое x1
        real fun = 0;
        real sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum += A[i * n + j] * X2[j];
                }
            }
            fun = F[i] - sum;
            X1[i] = fun / A[i * n + i];
            sum = 0;
        }
        return X1;

    }
}

void Relaxation_met(real *A, real *F, real *&X1, real *&X2, int n, int k, double w) {
    if (k % 2 == 1) {//нечетная итерация x1 прошлое x2
        real fun = 0;
        real sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j < i) {
                    sum += A[i * n + j] * X2[j];
                } else if (j > i) {
                    sum += A[i * n + j] * X1[j];
                }
            }
            fun = F[i] - sum;
            X2[i] = (1 - w) * X1[i] + w * fun / A[i * n + i];
            sum = 0;
        }
    } else {//четная итерация x2 прошлое x1
        real fun = 0;
        real sum = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (j < i) {
                    sum += A[i * n + j] * X1[j];
                } else if (j > i) {
                    sum += A[i * n + j] * X2[j];
                }
            }
            fun = F[i] - sum;
            X1[i] = (1 - w) * X1[i] + w * fun / A[i * n + i];
            sum = 0;
        }
    }
}

bool MaxCheck(real *X1, real *X2, int k, int n, double e) {
    if (k % 2 == 1) {//нечетная итерация x1 прошлое x2
        real maximum = fabs(X2[0] - X1[0]);
        for (int i = 1; i < n; i++) {
            maximum = max(maximum, fabs(X2[i] - X1[i]));
        }
        if (maximum < e) {
            return true;
        } else {
            return false;
        }
    } else {//четная итерация x2 прошлое x1
        real maximum = fabs(X1[0] - X2[0]);
        for (int i = 1; i < n; i++) {
            maximum = max(maximum, fabs(X1[i] - X2[i]));
        }
        if (maximum < e) {
            return true;
        } else {
            return false;
        }
    }
}

int main() {

    const double e = 0.0001;
    const int k_max = 1000;
    FillA(A_, n_);
    FillX(X_, n_, m_);
    cout << "X: ";
    for (int i = 0; i < 5; i++) {
        cout << X_[i] << ", ";
    }
    cout << endl << endl;

//    // Находим f=Ay
    MatrixMultiplication(F_, A_, n_, n_, X_, n_, 1);

    real *X0 = new real[n_];
    for (int i = 0; i < n_; i++) {
        X0[i] = F_[i] / A_[i * n_ + i];
    }


    //выполним метод Якоби
    int k = 1;
    bool new_iter = true;
    bool Is_k_max = false;
    real *X1 = new real[n_];
    real *X2 = new real[n_];
    for (int i = 0; i < n_; i++) {
        double x = X0[i];
        X1[i] = X0[i];
    }
    while (new_iter) {
        Jacobi_met(A_, F_, X1, X2, n_, k);
        while (MaxCheck(X1, X2, k, n_, e)) {
            new_iter = false;
            break;
        }
        while (k == k_max) {
            new_iter = false;
            Is_k_max = true;
            break;
        }
        k++;
    }
    k--;
    if (k % 2 == 1) {
        cout << "Jacobi_met: " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X2[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }


    } else {
        cout << "Jacobi_met: " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X1[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    }
    //закончили метод Якоби



    //выполним метод релаксации w=0,5
    k = 1;
    new_iter = true;
    Is_k_max = false;
    for (int i = 0; i < n_; i++) {
        X1[i] = X0[i];
    }
    while (new_iter) {
        Relaxation_met(A_, F_, X1, X2, n_, k, 0.5);
        while (MaxCheck(X1, X2, k, n_, e)) {
            new_iter = false;
            break;
        }
        while (k == k_max) {
            new_iter = false;
            Is_k_max = true;
            break;
        }
        k++;
    }
    k--;
    if (k % 2 == 1) {
        cout << "Relaxation_met (w=0,5): " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X2[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    } else {
        cout << "Relaxation_met (w=0,5): " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X1[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    }

    /////////////////////////////////////////////////////////////


    //выполним метод релаксации w=1
    k = 1;
    new_iter = true;
    Is_k_max = false;
    for (int i = 0; i < n_; i++) {
        X1[i] = X0[i];
    }
    while (new_iter) {
        Relaxation_met(A_, F_, X1, X2, n_, k, 1);
        while (MaxCheck(X1, X2, k, n_, e)) {
            new_iter = false;
            break;
        }
        while (k == k_max) {
            new_iter = false;
            Is_k_max = true;
            break;
        }
        k++;
    }
    k--;
    if (k % 2 == 1) {
        cout << "Relaxation_met (w=1): " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X2[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    } else {
        cout << "Relaxation_met (w=1): " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X1[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    }

    /////////////////////////////////////////////////////////



    //выполним метод релаксации w=1,5
    k = 1;
    new_iter = true;
    Is_k_max = false;
    for (int i = 0; i < n_; i++) {
        X1[i] = X0[i];
    }
    while (new_iter) {
        Relaxation_met(A_, F_, X1, X2, n_, k, 1.5);
        while (MaxCheck(X1, X2, k, n_, e)) {
            new_iter = false;
            break;
        }
        while (k == k_max) {
            new_iter = false;
            Is_k_max = true;
            break;
        }
        k++;
    }
    k--;
    if (k % 2 == 1) {
        cout << "Relaxation_met (w=1,5): " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X2[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    } else {
        cout << "Relaxation_met (w=1,5): " << endl;
        for (int i = 0; i < 5; i++) {
            cout << X1[i] << ", ";
        }
        if (Is_k_max) {
            cout << endl << "k = k_max: " << k << endl;
        } else {
            cout << endl << "k: " << k << endl;
        }
    }

    /////////////////////////////////////////////

}