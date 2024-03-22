//#pragma clang diagnostic push
//#pragma ide diagnostic ignored "cert-msc50-cpp"

#include <iostream>
#include <valarray>
#include <chrono>
#include <random>

using namespace std;
using real = double;
const int n = 100;
const int e=1;
const int k=1000;
real *A;
real * u;
real * v;
real lambda_1;
real lambda_2;

// Евклидова норма
double EuclideanNorm(real *X) {
    real sum = 0;
    for (int i = 0; i < n; i++) {
        sum += pow(X[i], 2);
    }
    return sqrt(sum);
}
  real MaxNorm(real *X){
    real max=abs(X[0]);
    for(int i=1;i<n;i++){
        if(abs(X[i])> max){
            max=abs(X[i]);
        }
    }
    return max;
}
real ScalarMultiplication(real *X,real *Y){
    real result=0;
    for(int i=0;i<n;i++){
        result+=X[i]*Y[i];
    }
    return result;
}

int FindMaxIndex(real* x){
    int res=0;
    for(int i=1;i<n;i++){
        if(abs(x[i])>abs(x[res])){
            res=i;
        }
    }
    return res;
}

void FillA() {
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
            A[0] = sum + pow(10, 2 - e);
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
void * MatrixMultiplication(real *&F, const real *A_,
                            int first_num_rows, int first_num_col,
                            const real *X,
                            int second_num_rows, int second_num_col) {
    for (int i = 0; i < first_num_rows; i++) {
        for (int j = 0; j < second_num_col; j++) {
            real elem = 0;
            for (int o = 0; o < first_num_col; o++) {
                elem += A_[i * first_num_col + o] * X[o * second_num_col + j];
            }
            F[i] = elem;
        }
    }
}
void Fill_U_zero(){
    u=new real[n];
    for(int i=0;i<n;i++){
        u[i]=1;
    }
}

void * Fill_U(real *&u_){
    double norm_v= MaxNorm(v);
    for(int i=0;i<n;i++){
        u[i]=v[i]/norm_v;
    }
}

void * Fill_V(){
    MatrixMultiplication(v,A,n,n,u,n,1);
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
void Find_lambda_sign(){
    int index= FindMaxIndex(v);
    lambda_1=v[index]*sign(u[index]);
}
void Find_lambda_scal(){
    lambda_2= ScalarMultiplication(v,u)/ ScalarMultiplication(u,u);
}
void Print_vec(real *x){
    cout<<"( ";
    for(int i=0;i<n;i++){
        cout<<x[i]<<", ";
    }
    cout<<")";
}


int main() {
    v=new real[n];
    u=new real[n];
    A=new real [n*n];
    FillA();
    Fill_U_zero();
    for(int i=0;i<k;i++){
        Fill_V();//перезаписываем v
        Find_lambda_sign();//перезаписываем lambda_1
        Find_lambda_scal();//перезаписываем lambda_2
        if(i!=k-1){
        Fill_U(u);//перезаписываем u, на последней итерации u не перезаписываем
        }
    }
    // находим вектор отклонения v-lambda*u
    real * vec_1=new real[n];
    real * vec_2=new real[n];
    for(int i=0;i<n;i++){
        vec_1[i]=v[i]-lambda_1*u[i];
        vec_2[i]=v[i]-lambda_2*u[i];
    }

    //Выводим результат

    cout<<"Result 1: with lambda located via sign\n";
    cout<<"u: ";
    Print_vec(u);
    cout<<"\nlambda_1 = "<<lambda_1;
    cout<<"\nEuclidean Norm : "<<EuclideanNorm(vec_1);
    cout<<"\nMax Norm : "<<MaxNorm(vec_1);

    cout<<"\nResult 2: with lambda found through the quotient of scalar products\n";
    cout<<"u: ";
    Print_vec(u);
    cout<<"\nlambda_2 = "<<lambda_2;
    cout<<"\nEuclidean Norm : "<<EuclideanNorm(vec_2);
    cout<<"\nMax Norm : "<<MaxNorm(vec_2);



}