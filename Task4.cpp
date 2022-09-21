
// #include "Task3.hpp"
#include <math.h>
#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"
#include "Creating_matrix.hpp"


int main(){

    // Define variables that wil be used in the jacobian eigenevalue solver
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    int maxiter = 10000;
    bool converged;

    int k;
    int l;
    double N = 6; // Matrix size

    arma::mat A = tri_matrix(N, false); // Define a matrix
    // Run the Jacobian eigenvalue solver
    jacobi_eigensolver(A, pow(10,-20), eigenvalues, eigenvectors, maxiter, iterations, converged);

    std::cout << eigenvalues; // prints eigenevalues
    std::cout << eigenvectors; // prints eigenvectors

    return 0;
}


