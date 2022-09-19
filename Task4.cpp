
// #include "Task3.hpp"
#include <math.h>
#include <armadillo>
#include <iostream>

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);


int main(){

    arma::mat A = arma::mat();
    A.load("fileT");
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    int maxiter;
    bool converged;

    int k;
    int l;

    // jacobi_eigensolver(A, pow(10,8), eigenvalues, eigenvectors, maxiter, iterations, converged);
    arma::mat R = arma::mat(A.n_rows, A.n_cols, arma::fill::eye);
    jacobi_rotate(A, R, k, l);

}


// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

    double tau = (A(l,l) - A(k,k))/(2*A(k,k));
    
    double tan = 1/(tau+sqrt(1+pow(tau,2)));

    if(tau<0){
        tan = -1/(-tau+sqrt(1+pow(tau,2)));
    }
    

    double cos = 1/sqrt(1+pow(tan,2));
    double sin = cos*tan;


    //updating elements
    arma::mat A_old = A;
    
    A(k,k) = A_old(k,k)*pow(cos,2) - 2*A_old(k,l)*cos*sin + 2*A_old(l,l)*pow(sin,2);
    A(l,l) = A_old(l,l)*pow(cos,2) + 2*A_old(k,l)*cos*sin + 2*A_old(k,k)*pow(sin,2);
    A(k,l) = 0;
    A(l,k) = 0;

    for(double i = 0; i < (A.n_cols-1); i++){
        if(i != k && i != l){
            A(i,k) = A_old(i,k)*cos - A_old(i,l)*sin;
            A(k,i) = A(i,k);
            A(i,l) = A_old(i,l)*cos - A_old(i,k)*sin;
            A(l,i) = A(i,l);
        }
    }

    arma::mat R_old = R;

    for(double i = 0; i < (R.n_cols-1); i++){
        R(i,k) = R_old(i,k)*cos - R(i,l)*sin;
        R(i,l) = R_old(i,l)*cos + R(i,k)*sin;
    }

    

}

// Jacobi method eigensolver:
// - Runs jacobo_rotate until max off-diagonal element < eps
// - Writes the eigenvalues as entries in the vector "eigenvalues"
// - Writes the eigenvectors as columns in the matrix "eigenvectors"
//   (The returned eigenvalues and eigenvectors are sorted using arma::sort_index)
// - Stops if it the number of iterations reaches "maxiter"
// - Writes the number of iterations to the integer "iterations"
// - Sets the bool reference "converged" to true if convergence was reached before hitting maxiter
void jacobi_eigensolver(arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    
    arma::mat R = arma::mat(A.n_rows, A.n_cols, arma::fill::eye);

    int k;
    int l;
    double max_value = max_offdiag_symmetric(A, k, l);
    //iterations = 1;

    // while(max_value > eps){
    //     jacobi_rotate(A, R, k, l);
    //     max_value = max_offdiag_symmetric(A, k, l);
    //     iterations += 1;
    // }

    //jacobi_rotate1(A, R, k, l);

    std::cout<< A;
    std::cout<< R;
    
    }

    double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){

    
    double max_value = 0;

    for(int i=0; i < A.n_rows-1; i++){
        for(int j=0;j < A.n_cols-1; j++){

            if(i!=j){
                if(abs(A(i,j))>max_value){
                 max_value = abs(A(i,j));
                 k=i;
                 l=j;
                }
            } 
        }
   }


    return max_value;
}