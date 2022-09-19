
// #include "Task3.hpp"
#include <math.h>
#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"
#include "Creating_matrix.hpp"


int main(){

    arma::mat A = tri_matrix(6, false);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    int maxiter;
    bool converged;

    int k;
    int l;

    jacobi_eigensolver(A, pow(10,-20), eigenvalues, eigenvectors, maxiter, iterations, converged);

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

    return 0;
}


