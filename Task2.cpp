#include <armadillo>
#include <iostream>
#include <math.h>

int main(){
 
    int N = 6;
    int n = N;
    int m = N;
    int num_steps = N+1;
    double h = 1/num_steps;   //Our step_size
    double a = -1/pow(h,2);
    double b = 2/pow(h,2);


    arma::mat A = arma::mat(n,m).fill(0);
    b = A(0,0);
    a = A(0,1);
    b = A(N-1, N-1);
    a=A(N-1, N-2);
    for(int i=1; i < N-2; i++){
        for(int j=1;j<N-2; j++){

            if(i==j){
                b=A(i,j);
                a=A(i-1,j);
                a=A(i+1,j);
            }

        }

    }

arma::vec eigenval;
arma::mat eigenvec;

//arma::eig_sym(eigenval, eigenvec, A);
//std::cout << eigenvec;
//std::cout << eigenval;
return 0;
     

}