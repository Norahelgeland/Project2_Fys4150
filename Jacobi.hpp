#include <math.h>
#include <armadillo>
#include <iostream>
//#include "Max_offdiag_symmetric.hpp"

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, arma::mat& eigenvectors, 
                        const int maxiter, int& iterations, bool& converged);
double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);


void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

    double tau = (A(l,l) - A(k,k))/(2*A(k,l)); //Define tau

    // Define tan with dependence on the value of tau
    double tan = 1/(tau+sqrt(1+pow(tau,2)));
    if(tau<0){
        tan = -1/(-tau+sqrt(1+pow(tau,2)));
    }
    // Define cos and sin
    double cos = 1/sqrt(1+pow(tan,2));
    double sin = cos*tan;

    //Now we are starting to update the elements

    arma::mat A_old = A; // Save the matrix before updating

    A(k,k) = A_old(k,k)*pow(cos,2) - 2*A_old(k,l)*cos*sin + A_old(l,l)*pow(sin,2);
    A(l,l) = A_old(l,l)*pow(cos,2) + 2*A_old(k,l)*cos*sin + A_old(k,k)*pow(sin,2);
    A(k,l) = 0;
    A(l,k) = 0;

    // Iterate trhoigh every i, except i = k, and l and update A.
    for(double i = 0; i < (A.n_cols); i++){
        if(i != k && i != l){
            A(i,k) = A_old(i,k)*cos - A_old(i,l)*sin;
            A(k,i) = A(i,k);
            A(i,l) = A_old(i,l)*cos - A_old(i,k)*sin;
            A(l,i) = A(i,l);
        }
    }
    // Update the matrix R
    arma::mat R_old = R; // Save the matrix before updating

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
    arma::mat D = A; // Define a new matrix that wil become the diagonal matrix
    int k;
    int l;
    double max_value = max_offdiag_symmetric(D, k, l); // Get the max value of D
    iterations = 0;

    // Run the jacobian rotation function untill the criteria is met
    while(max_value > pow(10,-10)){
        jacobi_rotate(D, R, k, l);
        max_value = max_offdiag_symmetric(D, k, l);
        iterations += 1;
    }
    // Pick out the eigenevalues and eigenevectors
    eigenvalues = D.diag();
    eigenvectors = R;

    }

    // Function that returns the bigest absolute value (that is not on the diagonal) 
    // and changes k and l to be the position of said element
    double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){

    double max_value = 0; // Define a variable that wil hold the max value
 
    // Runs trough every position in the matrix and find the maximal value
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
