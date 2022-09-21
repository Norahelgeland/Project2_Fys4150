#include <armadillo>
#include <iostream>
#include <math.h>
#include "Creating_matrix.hpp"

int main(){
    
    // Define dimentions of matrix
    double N = 6;
    // Define matrix
    arma::mat A = tri_matrix(N, false);

    A.save("file1"); // Save the matrix
    // Define vector and matrix that wil hold the eigenvectors and eigenvalues
    arma::vec eigenval;
    arma::mat eigenvec;

    arma::eig_sym(eigenval, eigenvec, A); // Find eigenevalues and eigenvectors
    std::cout << eigenvec;
    std::cout << eigenval;
    arma::mat V = arma::normalise(eigenvec); // Normalize the eigenvectors
    std::cout << V;

    return 0;     

}