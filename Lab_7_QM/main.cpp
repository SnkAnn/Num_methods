//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "cert-msc50-cpp"

#include <iostream>
#include <valarray>
#include <chrono>
#include <random>

using namespace std;
using real = float;
const int n = 4;
real *A;
real *M;
real *min_M;

float randomFloat(int a, int b)
{
    return (float)(a + (rand() % (b - a))) + (float)(rand()) / (float)(RAND_MAX);
}

void FillA() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] =randomFloat(-50,50);
        }
    }
}
real * MatrixMultiplication( const real *A_,
                             int first_num_rows, int first_num_col,
                             const real *X,
                             int second_num_rows, int second_num_col) {

    real * F = new real[first_num_rows*second_num_col];
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                elem += A_[i * first_num_col + k] * X[k * second_num_col + j];
            }
            F[i*n+j] = elem;
        }
    }
    return F;
}

void * MatrixMultiplication(real *&F, const real *A_,
                            int first_num_rows, int first_num_col,
                            const real *X,
                            int second_num_rows, int second_num_col) {
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int k = 0; k < first_num_col; k++) {
                elem += A_[i * first_num_col + k] * X[k * second_num_col + j];
            }
            F[i*n+j] = elem;
        }
    }
}
void Write_sqr_Matrix(int size,real *Mat){
    for(int i=0;i<size;i++){
        cout<<"[ ";
        for(int j=0;j<size;j++){
            cout<<Mat[size*i+j]<<" ";
        }
        cout<<"]"<<endl;
    }
}
void Fill_M_k_part(int k,real *&M){
    for(int i=0;i<n;i++){
        if(i!=k){
            for(int j=0;j<n;j++){
                if(i!=j){
                    M[k*n*n+i*n+j]=0;
                }else{
                    M[k*n*n+i*n+j]=1;
                }
            }

        }else{
            for(int j=0;j<n;j++){
                if(j!=k){
                M[k*n*n+i*n+j]=-A[n*(k+1)+j]/A[n*(k+1)+k];
                }else{
                    M[k*n*n+i*n+j]=1/A[n*(k+1)+k];
                }
            }
        }
    }
}
void Fill_min_M_k_part(int k,real *&min_M){
    for(int i=0;i<n;i++){
        if(i!=k){
            for(int j=0;j<n;j++){
                if(i!=j){
                    min_M[k*n*n+i*n+j]=0;
                }else{
                    min_M[k*n*n+i*n+j]=1;
                }
            }
        }else{
            for(int j=0;j<n;j++){
                min_M[k*n*n+i*n+j]=A[n*(k+1)+j];
            }
        }
    }
}
int main() {
    srand(time(0));
    A=new real [n*n];
    FillA();
    //ищем след
    real sum=0;
    for(int i=0;i<n;i++){
        sum+=A[i*n+i];
    }
    cout<<"Matrix A:"<<endl;
    Write_sqr_Matrix(n,A);
    M=new real [n*n*(n-1)];
    min_M=new real[n*n*(n-1)];
    int k=n-2;
    real a;
    real *min_M_k=new real[n*n];
    real *M_k=new real[n*n];
    while(k>=0){
        a=A[n*(k+1)+k];
        if(a!=0 && abs(a)>=pow(10,-8)){
            Fill_M_k_part(k,M);// заполняем матрицу M
            Fill_min_M_k_part(k,min_M);// заполняем матрицу M^-1
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    min_M_k[i*n+j]=min_M[k*n*n+n*i+j];
                    M_k[i*n+j]=M[k*n*n+n*i+j];
                }
            }
//            cout<<endl;
//            Write_sqr_Matrix(n,min_M_k);
//            cout<<endl;
//            Write_sqr_Matrix(n,M_k);
            MatrixMultiplication(A,min_M_k,n,n,
                                 MatrixMultiplication(A,n,n,M_k,n,n),
                                 n,n);//находим матрицу A на k-ом шаге

//            cout<<k+"   "<<endl;
//            Write_sqr_Matrix(n,A);
        }else{
            cout<<"The leading element turned out to be too small, let's rewrite the matrix A."<<endl;
            FillA();
            cout<<"Matrix A:"<<endl;
            Write_sqr_Matrix(n,A);
            //ищем след
            sum=0;
            for(int i=0;i<n;i++){
                sum+=A[i*n+i];
            }
            k=n-2;
        }
        k--;
    }
    //все матрицы M хранятся в массиве матриц M под номерами в массиве, пример M[n-1] храниться на месте n-2, а M[1] на месте 0

    cout<<endl<<"Frobenius matrix:"<<endl;
    Write_sqr_Matrix(n,A);
    cout<<endl<<"M matrix:"<<endl;
    for(int i=n-1;i>0;i--){
        cout<<"M"<<i<<" :"<<endl;
        for(int j=0;j<n;j++){
            cout<<"[ ";
            for(int t=0;t<n;t++){
                cout<<M[(i-1)*n*n+j*n+t] <<" ";
            }
            cout<<"]"<<endl;
        }
        cout<<endl;
    }
    cout<<"p1: "<<A[0]<<endl;
    cout<<"Sp A : "<<sum;
}