#include <armadillo>
#include <iostream>
#include <math.h>

int main(){
 
    double N = 6;
    double n = N;
    double m = N;
    double num_steps = N;
    double h = 1/num_steps;   //Our step_size
    double a = -1/pow(h,2);
    double b = 2/pow(h,2);
    arma::mat A = arma::mat(n,m).fill(0.);
 
    A(0,0)=b;
    A(0,1)=a;
    A(N-1, N-1)=b;
    A(N-1, N-2)=a;

    for(int i=1; i < N-1; i++){
        for(int j=1;j<N-1; j++){

            if(i==j){
                A(i,j)=b;
                A(i,j-1)=a;
                A(i,j+1)=a;
                
            }      
        }
   }

A.save("file1");
arma::vec eigenval;
arma::mat eigenvec;

arma::eig_sym(eigenval, eigenvec, A);
std::cout << eigenvec;
std::cout << eigenval;
std::cout << A;
arma::mat V = arma::normalise(eigenvec); 
std::cout << V;
return 0;     

}