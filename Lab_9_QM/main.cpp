//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "cert-msc50-cpp"

#include <iostream>
#include <valarray>
#include <chrono>
#include <random>

using namespace std;
using real = double;
const int n = 4;
real *A;
real *M;
real *min_M;
real *P;
real x0=-50;
real x1=50;
real xk;
real road;
const real e=0.00001;

float randomFloat(int a, int b) {
    return (float) (a + (rand() % (b - a))) + (float) (rand()) / (float) (RAND_MAX);
}

void FillA() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i*n+j]= (-50 + (rand() % (50 + 50))) +  (rand()) / (RAND_MAX);
           // A[i * n + j] = randomFloat(-50, 50);
        }
    }
}

real *MatrixMultiplication(const real *A_,
                           int first_num_rows, int first_num_col,
                           const real *X,
                           int second_num_rows, int second_num_col) {

    real *F = new real[first_num_rows * second_num_col];
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                elem += A_[i * first_num_col + k] * X[k * second_num_col + j];
            }
            F[i * n + j] = elem;
        }
    }
    return F;
}

void *MatrixMultiplication(real *&F, const real *A_,
                           int first_num_rows, int first_num_col,
                           const real *X,
                           int second_num_rows, int second_num_col) {
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                elem += A_[i * first_num_col + k] * X[k * second_num_col + j];
            }
            F[i * n + j] = elem;
        }
    }
}

void Write_sqr_Matrix(int size, real *Mat) {
    for (int i = 0; i < size; i++) {
        cout << "[ ";
        for (int j = 0; j < size; j++) {
            cout << Mat[size * i + j] << " ";
        }
        cout << "]" << endl;
    }
}

void Fill_M_k_part(int k, real *&M) {
    for (int i = 0; i < n; i++) {
        if (i != k) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    M[k * n * n + i * n + j] = 0;
                } else {
                    M[k * n * n + i * n + j] = 1;
                }
            }

        } else {
            for (int j = 0; j < n; j++) {
                if (j != k) {
                    M[k * n * n + i * n + j] = -A[n * (k + 1) + j] / A[n * (k + 1) + k];
                } else {
                    M[k * n * n + i * n + j] = 1 / A[n * (k + 1) + k];
                }
            }
        }
    }
}

void Fill_min_M_k_part(int k, real *&min_M) {
    for (int i = 0; i < n; i++) {
        if (i != k) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    min_M[k * n * n + i * n + j] = 0;
                } else {
                    min_M[k * n * n + i * n + j] = 1;
                }
            }
        } else {
            for (int j = 0; j < n; j++) {
                min_M[k * n * n + i * n + j] = A[n * (k + 1) + j];
            }
        }
    }
}

void FindP(real *&X) {
   // srand(time(0));
    A = new real[n * n];
    FillA();
    //ищем след
    real sum = 0;
    for (int i = 0; i < n; i++) {
        sum += A[i * n + i];
    }
    M = new real[n * n * (n - 1)];
    min_M = new real[n * n * (n - 1)];
    int k = n - 2;
    real a;
    real *min_M_k = new real[n * n];
    real *M_k = new real[n * n];
    while (k >= 0) {
        a = A[n * (k + 1) + k];
        if (a != 0 && abs(a) >= pow(10, -8)) {
            Fill_M_k_part(k, M);// заполняем матрицу M
            Fill_min_M_k_part(k, min_M);// заполняем матрицу M^-1
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    min_M_k[i * n + j] = min_M[k * n * n + n * i + j];
                    M_k[i * n + j] = M[k * n * n + n * i + j];
                }
            }
            MatrixMultiplication(A, min_M_k, n, n,
                                 MatrixMultiplication(A, n, n, M_k, n, n),
                                 n, n);//находим матрицу A на k-ом шаге

        } else {
            cout << "The leading element turned out to be too small, let's rewrite the matrix A." << endl;
            FillA();
            //ищем след
            sum = 0;
            for (int i = 0; i < n; i++) {
                sum += A[i * n + i];
            }
            k = n - 2;
        }
        k--;
    }
    //все матрицы M хранятся в массиве матриц M под номерами в массиве, пример M[n-1] храниться на месте n-2, а M[1] на месте 0
    for (int i = 0; i < n; i++) {
        X[i] = A[i];
    }

}
void FindDerivative(real *X,int sizeX,  real *&Y,int sizeY){
    int s=sizeX;
    int l=sizeY;
    int step=1;
    for(int i=s-1;i>0;i--){
        real x=X[i-1];
        Y[l-step]=X[i-1]*step;
        step++;
    }
    if(s==l){
    Y[0]=step;
    }
}
real  FindCountDer(real x,real *X){
    real sum=pow(x,n-1);
    for(int i=1;i<n-1;i++){
        sum+=pow(x,n-i)*X[i-1];
    }
    sum+=X[n-1];
    return sum;
}
int sign(real x){
    if(x==0){
        return 0;
    }else if(x>0){
        return 1;
    }else{
        return -1;
    }
}
real FindF(real xk, int size, real *X){
    int s=size;
    real sum=0;
    int step=1;
    sum+=X[s-1];
    for(int i=s-2;i>=0;i--){
        sum+=X[i]*pow(xk,step);
        step++;
    }
    return sum;
}
void MedianMethod(real &xk, real &delta, real &road,real *X){
    delta/=2;
    real lastx=xk;
xk-=delta* sign(FindF(xk,n,X));
road= abs(xk-lastx);
}
void NewtonMethod(real &xk, real* F, real *F_der){
    xk-= (FindF(xk, n, F) / FindF(xk, n - 1, F_der));
}
void ConstDerMethod(real &xk, real* F, real F_der_zero){
    xk-= (FindF(xk, n, F) /  F_der_zero);
}
void SecantsMethod(real &xk,real &xkLater,real*F){
    real laterXK=xk;
    xk-=FindF(xk, n, F)*((xk-xkLater)/(FindF(xk, n, F)-FindF(xkLater, n, F)));
    xkLater=laterXK;
}
int main() {
    P = new real[n];
    FindP(P);
    for(int i=0;i<n;i++){
        P[i]=-P[i];
    }
    cout<<"P: ";
    for(int i=0;i<n;i++){
        cout<<" "<<P[i];
    }
    cout<<"\n";

    real*P_one=new real[n];
    real *P_second=new real[n-1];
    real *Der_P_first;
    FindDerivative(P,n,P_one,n);
   // P_one=P;
    FindDerivative(P_one,n,P_second,n-1);

    cout<<"P_one: ";
    for(int i=0;i<n;i++){
        cout<<" "<<P_one[i];
    }
    cout<<"\n";

    cout<<"P_second: ";
    for(int i=0;i<n;i++){
        cout<<" "<<P_second[i];
    }
    cout<<"\n";
    //все для 2
    real delta=abs(x1-x0)/2;
    xk=(x0+x1)/2;
    road=abs(x1-x0);
    int k=2;
    while(road>2*e){
        MedianMethod(xk,delta,road,P_one);
        k++;
    }
    cout<<"\ncount of iters = "<<k<<"\n x median = "<<xk;
    xk=x0;
    k=1;
    while(abs(FindF(xk,n,P_one))>e){
        NewtonMethod(xk,P_one,P_second);
        int xe=xk;
        k++;
    }
    cout<<"\ncount of iters = "<<k<<"\n x newton = "<<xk;

    k=1;
    xk=x0;
    real ConstDer=FindF(xk, n - 1, P_second);
    while(abs(FindF(xk,n,P_one))>e){
        ConstDerMethod(xk,P_one,ConstDer);
        int xe=xk;
        k++;
    }
    cout<<"\ncount of iters = "<<k<<"\n x constDer = "<<xk;

    k=2;
    real xkLater=x0;
    xk=x1;
    while(abs(FindF(xk,n,P_one))>e){
        SecantsMethod(xk,xkLater,P_one);
        int xe=xk;
        k++;
    }
    cout<<"\ncount of iters = "<<k<<"\n x secant = "<<xk;


    //часть 2

    cout<<"\n";

    real* y=new real[n];
    real lambda=xk;
    int t=0;
    for(int i=n-1;i>=0;i--){
        y[t]=pow(lambda,i);
        t++;
    }
    real *S=new real[n*n];
    for(int i=0;i<n;i++){
         S[i*n+i]=1;
    }
    for(int i=n-2;i>=0;i--){
        real *M_part=new real[n*n];
        for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
                M_part[j*n+k]=M[i*n*n+j*n+k];
            }
        }
        MatrixMultiplication(S,S,n,n,M_part,n,n);
    }
    real *u=new real[n];
    MatrixMultiplication(u,S,n,n,y,n,1);
    real *numMult=new real[n];
    real *Ay=new real[n];
    MatrixMultiplication(Ay,A,n,n,u,n,1);
    for(int i=0;i<n;i++){
        numMult[i]=numMult[i]-lambda*u[i];
    }
    cout<<"Cheking Au-lu=0: ";
    for(int i=0;i<n;i++){
        cout<<" "<<numMult[i];
    }
    cout<<"\n";
}