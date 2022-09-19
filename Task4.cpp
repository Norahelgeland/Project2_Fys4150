
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
    A.load("file4");
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    int maxiter;
    bool converged;

    int k;
    int l;
    //max_offdiag_symmetric(A, k, l);
    jacobi_eigensolver(A, pow(10,-20), eigenvalues, eigenvectors, maxiter, iterations, converged);
    //arma::mat R = arma::mat(A.n_rows, A.n_cols, arma::fill::eye);


}


// Performs a single Jacobi rotation, to "rotate away"
// the off-diagonal element at A(k,l).
// - Assumes symmetric matrix, so we only consider k < l
// - Modifies the input matrices A and R
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

    double tau = (A(l,l) - A(k,k))/(2*A(k,l));
    
    double tan = 1/(tau+sqrt(1+pow(tau,2)));

    if(tau<0){
        tan = -1/(-tau+sqrt(1+pow(tau,2)));
    }
    

    double cos = 1/sqrt(1+pow(tan,2));
    double sin = cos*tan;


    //updating elements
    arma::mat A_old = A;
    
    A(k,k) = A_old(k,k)*pow(cos,2) - 2*A_old(k,l)*cos*sin + A_old(l,l)*pow(sin,2);
    A(l,l) = A_old(l,l)*pow(cos,2) + 2*A_old(k,l)*cos*sin + A_old(k,k)*pow(sin,2);
    A(k,l) = 0;
    A(l,k) = 0;

    for(double i = 0; i < (A.n_cols); i++){
        if(i != k && i != l){
            A(i,k) = A_old(i,k)*cos - A_old(i,l)*sin;
            A(k,i) = A(i,k);
            A(i,l) = A_old(i,l)*cos - A_old(i,k)*sin;
            A(l,i) = A(i,l);
        }
    }

    arma::mat R_old = R;

    for(double i = 0; i < (R.n_cols); i++){
        R(i,k) = R_old(i,k)*cos - R_old(i,l)*sin;
        R(i,l) = R_old(i,l)*cos + R_old(i,k)*sin;
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
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged){
    
    arma::mat R = arma::mat(A.n_rows, A.n_cols, arma::fill::eye);
    arma::mat D = A;
    int k;
    int l;
    double max_value = max_offdiag_symmetric(D, k, l);
    iterations = 0;

    //std::cout<< max_value;
    //std::cout<< eps;

    while(max_value > pow(10,-8)){
        jacobi_rotate(D, R, k, l);
        max_value = max_offdiag_symmetric(D, k, l);
        std::cout<< l;
        std::cout<< "#################################################################################################################################################";
        std::cout<< k;
        iterations += 1;
    }

    //std::cout<< iterations; 
    std::cout<< D;
    // //std::cout<< A;
    std::cout<< "#################################################################################################################################################";
    std::cout<< R;

    std::cout<< "#################################################################################################################################################";

    arma::vec eigenvaltrue;
    arma::mat eigenvectrue;
    arma::eig_sym(eigenvaltrue, eigenvectrue, A);
    std::cout << eigenvectrue;
    std::cout << eigenvaltrue;

    std::cout<< "#################################################################################";
    std::cout<< A;


    }

    double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){

    
    double max_value = 0;
 
    for(int i=0; i < A.n_rows; i++){
        for(int j=0;j < A.n_cols; j++){

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