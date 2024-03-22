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
const double e=0.0000001;//погрешность
real *A_;// исходная матрица
real *X_;//точное решение
real *F_;

void FillA(real *&A, int n, int m, int k) {
    A = new real[n*n];
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

void FillX(real *&X, int n, int m) {
    X = new real[n];
    for (int i = 0; i < n; i++) {
        X[i] = m;
        m++;
    }
}

void * MatrixMultiplication(real *&F, const real *A,
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
real * MatrixMultiplication( const real *A,
                            int first_num_rows, int first_num_col,
                            const real *X,
                            int second_num_rows, int second_num_col) {

   real * F = new real[first_num_rows];
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                elem += A[i * first_num_col + k] * X[k * second_num_col + j];
            }
            F[i] = elem;
        }
    }
    return F;
}
real ScalarMultiplication(real *X,real *Y,int n){
    real result=0;
    for(int i=0;i<n;i++){
        result+=X[i]*Y[i];
    }
    return result;
}
real Find_alf(real scal_R_gr,real*P_gr,real *A_p){
    return scal_R_gr / ScalarMultiplication(A_p,P_gr,n_);
}
real Find_bet(real *R_gr,real scal_R_gr ){
    return ScalarMultiplication(R_gr,R_gr,n_)/scal_R_gr;
}
void Find_X(real *&X,real alf,real*P){
    for(int i=0;i<n_;i++){
        X[i]+=alf*P[i];
        double x=X[i];
    }
}
void Find_R(real *&R,real alf,real*A_p){
    for(int i=0;i<n_;i++){
       R[i]-=alf*A_p[i];

    }
}
void Find_P(real*&P,real *R,real bet){
    for(int i=0;i<n_;i++){
        P[i]=R[i]+bet*P[i];
    }
}
void Grad_met(real *&X_gr,real*&R_gr,real *&P_gr,real &alf,real &bet,real &scal_R,real *&A_p){
    Find_X(X_gr,alf,P_gr);
    Find_R(R_gr,alf,A_p);
    bet=Find_bet(R_gr,scal_R);
    Find_P(P_gr,R_gr,bet);
    A_p=MatrixMultiplication(A_,n_,n_,P_gr,n_,1);
    scal_R= ScalarMultiplication(R_gr,R_gr,n_);
    alf= Find_alf(scal_R,P_gr,A_p);
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
real RelativeError( real *X, const real *X_calculated, int n) {
    real relativeError;
    real *Distinction = new real[n];
    for (int i = 0; i < n; i++) {
        Distinction[i] = X[i] - X_calculated[i];
    }
    relativeError = EuclideanNorm(Distinction, n) / EuclideanNorm(X, n);
    return relativeError;
}
//считает невязку
double Residual(real *F,real *X_new){
    real* A_X= MatrixMultiplication(A_,n_,n_,X_new,n_,1);
    real* Res=new real[n_];
    for(int i=0;i<n_;i++){
        Res[i]=F[i]-A_X[i];
    }
    return EuclideanNorm(Res,n_);
}
void A_to_LDLt(real *&A, int n) {
    real *D = new real[n];
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            D[j] = A[j * n + i];
            A[j * n + i] /= A[i * n + i];
            for (int k = i + 1; k <= j; ++k) {
                A[j * n + k] -= A[j * n + i] * D[k];
            }
        }
    }
}
real *FillY(real *&Y, real *A, real *B, int n) {
    Y = new real[n];
    Y[0] = B[0];
    for (int i = 1; i < n; i++) {
        real sum = B[i];
        for (int k = 0; k <= i - 1; k++) {
            if (i > k) {
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
        Z[i] = Y[i] / A[i * n + i];
    }
    return Z;
}
real *FillX_find(real *&X, real *A, real *Z, int n) {
    X = new real[n];
    X[n - 1] = Z[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        double z = Z[i];
        real sum = z;
        for (int k = i + 1; k < n; k++) {
            if (k > i) {
                sum -= X[k] * A[k * n + i];
            }
        }
        X[i] = sum;
    }
}
int main() {
    const int l_max=50;
    int l=1;
    FillA(A_, n_, m_, k_);
    FillX(X_, n_,m_);
    // Находим f=Ay
    MatrixMultiplication(F_, A_, n_, n_, X_, n_, 1);

    //Все что надо для метода сод.градиентов
    real *X_gr=new real[n_];
    real *R_gr=new real[n_];
    real *P_gr=new real[n_];
    real alf;
    real bet;

    auto start_GR_method = chrono::high_resolution_clock::now();

    //пропишем все для нулевой итерации
    for(int i=0;i<n_;i++){
        R_gr[i]=F_[i];//так как вектор X при нулевой итерации нулевой
        P_gr[i]=R_gr[i];
    }
    real scal_R= ScalarMultiplication(R_gr,R_gr,n_);
    real* A_p=new real[n_];
    A_p=MatrixMultiplication(A_,n_,n_,P_gr,n_,1);
    alf= Find_alf(scal_R,P_gr,A_p);
    while(l<=l_max && EuclideanNorm(R_gr,n_)>=e){
        Grad_met(X_gr,R_gr,P_gr,alf,bet,scal_R,A_p);
        l++;
    }
    auto end_GR_method = chrono::high_resolution_clock::now();
    int time_GR_method = chrono::duration_cast<chrono::milliseconds>(end_GR_method - start_GR_method).count();
    real res_GR=Residual(F_,X_gr);//посчитали невязку
    real rel_error_GR= RelativeError(X_,X_gr,n_);

    //вычислили все по методу градиента

    //выполним LDLt метод
    // Решаем по методу LDLt_____________________________________

    real *Y = new real[n_];
    real *Z = new real[n_];
    real *Slay_X = new real[n_];// полученный X в LDLt
    //перепишем старую матрицу A, чтоб норма была нормальной
    real *A_old=new real[n_*n_];
    for(int i=0;i<n_;i++){
        for(int j=0;j<n_;j++){
            A_old[i*n_+j]=A_[i*n_+j];
        }
    }
    //расставили границы для времени
    auto start_LDLt = chrono::high_resolution_clock::now();
    A_to_LDLt(A_, n_);
    FillY(Y, A_, F_, n_);
    FillZ(Z, Y, A_, n_);
    FillX_find(Slay_X, A_, Z, n_);
    auto end_LDLt = chrono::high_resolution_clock::now();
    int time_LDLt = chrono::duration_cast<chrono::milliseconds>(end_LDLt - start_LDLt).count();

    //вернули старое A из записанного A_old для невязки
    for(int i=0;i<n_;i++){
        for(int j=0;j<n_;j++){
            A_[i*n_+j]=A_old[i*n_+j];
        }
    }
    real res_LDLt=Residual(F_,Slay_X);//посчитали невязку
    real rel_error_LDLt= RelativeError(X_,Slay_X,n_);
    // Конец решения по методу LDLt_____________________________


    cout << "X: ";
    for (int i = 0; i < 5; i++) {
        cout << X_[i] << ", ";
    }
    cout << endl << endl;
    cout<<"Conjugate gradient method:"<<endl;
    cout<<"X_gr:";
    for (int i = 0; i < 5; i++) {
        cout << X_gr[i] << ", ";
    }
    cout<<endl<<endl<<"l: "<<l--<<endl;
    cout<<"Norma: "<<res_GR<<endl<<"Relative error: "<<rel_error_GR<<endl<<"Time: "<<time_GR_method<<endl<<endl;

    cout<<" LDLt method:"<<endl;
    cout<<"X_LDLt:";
    for (int i = 0; i < 5; i++) {
        cout << Slay_X[i] << ", ";
    }
    cout<<endl;
    cout<<"Norma: "<<res_LDLt<<endl<<"Relative error: "<<rel_error_LDLt<<endl<<"Time: "<<time_LDLt<<endl;

}
// Вывод ответа в формате:
//1. Первые 5 координат вектора приближённого решения x*.
//невязка
//2. Относительная погрешность
//3. Время выполнения (можно приблизительно).