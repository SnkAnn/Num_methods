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
real *M_=new real[n_*n_];
real *M_min_one=new real[n_*n_];
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
                sum -= A[i * n + j];
                j++;
            }
            A[0] = sum + pow(10, 2 - k);
        } else {
            real sum = 0;
            int j = 0;
            while (j != n) {
                if (j != i) {
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
real Find_bet_pred(real *R_gr,real scal_R_gr){
    return ScalarMultiplication(MatrixMultiplication(M_min_one,n_,n_,R_gr,n_,1),R_gr,n_)/scal_R_gr;
}
void Find_X(real *&X,real alf,real*P){
    for(int i=0;i<n_;i++){
        X[i]+=alf*P[i];
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
void Find_P_pred(real *&P,real *R,real bet){
    real *Mr=new real[n_*n_];
    MatrixMultiplication(Mr,M_min_one,n_,n_,R,n_,1);
    for(int i=0;i<n_;i++){
        P[i]=Mr[i]+bet*P[i];
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
void Grad_met_pred(real *&X_gr,real*&R_gr,real *&P_gr,real &alf,real &bet,real &scal_R,real *&A_p){
    Find_X(X_gr,alf,P_gr);
    Find_R(R_gr,alf,A_p);
    bet=Find_bet_pred(R_gr,scal_R);
    Find_P_pred(P_gr,R_gr,bet);
    A_p=MatrixMultiplication(A_,n_,n_,P_gr,n_,1);
    scal_R= ScalarMultiplication(MatrixMultiplication(M_min_one,n_,n_,R_gr,n_,1),R_gr,n_);
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
void Fill_d(real *&d){
    real sum;
    for(int i=0;i<n_;i++){
        sum=0;
        for(int j=0;j<n_;j++){
            sum+=pow(A_[i*n_+j],2);
        }
        d[i]=sqrt(sum);
    }
}
void Fill_M_to_a(real *&M){
    for(int i=0;i<n_;i++){
        M[i*n_+i]=A_[i*n_+i];
        real a=A_[i*n_+i];
    }
}
void Fill_M_to_d(real *&M,real *d){
    for(int i=0;i<n_;i++){
        M[i*n_+i]=d[i];
    }
}
void Fill_M_minus(real *&M_min){
    for(int i=0;i<n_;i++){
        M_min[i*n_+i]=1/M_[i*n_+i];
    }
}
void pred_Grad_metod(real *&X,int l_max,int &l){
    real *R=new real[n_];
    real *P=new real[n_];
    real alf;
    real bet;
    //пропишем все для нулевой итерации
    for(int i=0;i<n_;i++){
        R[i]=F_[i];//так как вектор X при нулевой итерации нулево
    }
    MatrixMultiplication(P,M_min_one,n_,n_,R,n_,1);//находим P0
    real scal_R_m= ScalarMultiplication(MatrixMultiplication(M_min_one,n_,n_,R,n_,1),R,n_);
    real* A_p=new real[n_];
    A_p=MatrixMultiplication(A_,n_,n_,P,n_,1);
    alf= Find_alf(scal_R_m,P,A_p);
    while(l<=l_max && EuclideanNorm(R,n_)>=e){
        Grad_met_pred(X,R,P,alf,bet,scal_R_m,A_p);
        l++;
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


    //предобусловленный метод сопряженных градиентов
    //1) для M по a[i,i]
    real *X_gr_a=new real[n_];
    int l_a=1;
    auto start_predGR_a = chrono::high_resolution_clock::now();
    Fill_M_to_a(M_);
    Fill_M_minus(M_min_one);
    pred_Grad_metod(X_gr_a,l_max,l_a);
    auto end_predGR_a= chrono::high_resolution_clock::now();
    int time_predGR_a = chrono::duration_cast<chrono::milliseconds>(end_predGR_a - start_predGR_a).count();
    real res_predGR_a=Residual(F_,X_gr_a);//посчитали невязку
    real rel_error_predGR_a= RelativeError(X_,X_gr_a,n_);

    //2) для M по d[i,i]
    int l_d=1;
    real *d=new real[n_];
    real *X_gr_d=new real[n_];
    auto start_predGR_d = chrono::high_resolution_clock::now();
    Fill_d(d);
    Fill_M_to_d(M_,d);
    Fill_M_minus(M_min_one);
    pred_Grad_metod(X_gr_d,l_max,l_d);
    auto end_predGR_d= chrono::high_resolution_clock::now();
    int time_predGR_d = chrono::duration_cast<chrono::milliseconds>(end_predGR_d - start_predGR_d).count();
    real res_predGR_d=Residual(F_,X_gr_d);//посчитали невязку
    real rel_error_predGR_d= RelativeError(X_,X_gr_d,n_);
    //конец предобусловленного метода

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
    cout<<endl<<"l: "<<l--<<endl;
    cout<<"Norma: "<<res_GR<<endl<<"Relative error: "<<rel_error_GR<<endl<<"Time: "<<time_GR_method<<endl<<endl;

    ////////////////
    cout<<"Conjugate predefined gradient method with a[i*i]:"<<endl;
    cout<<"X_gr_a:";
    for (int i = 0; i < 5; i++) {
        cout << X_gr_a[i] << ", ";
    }
    cout<<endl<<"l: "<<l_a--<<endl;
    cout<<"Norma: "<<res_predGR_a<<endl<<"Relative error: "<<rel_error_predGR_a<<endl<<"Time: "<<time_predGR_a<<endl<<endl;

    //////////////////
    cout<<"Conjugate predefined gradient method with d[i*i]:"<<endl;
    cout<<"X_gr_d:";
    for (int i = 0; i < 5; i++) {
        cout << X_gr_d[i] << ", ";
    }
    cout<<endl<<"l: "<<l_d--<<endl;
    cout<<"Norma: "<<res_predGR_d<<endl<<"Relative error: "<<rel_error_predGR_d<<endl<<"Time: "<<time_predGR_d<<endl<<endl;


}
// Вывод ответа в формате:
//1. Первые 5 координат вектора приближённого решения x*.
//невязка
//2. Относительная погрешность
//3. Время выполнения (можно приблизительно).