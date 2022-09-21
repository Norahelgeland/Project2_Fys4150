
#include <math.h>
#include <armadillo>
#include <iostream>
#include "Jacobi.hpp"
#include "Creating_matrix.hpp"


int main(){


    double num_steps = 100;
    double N = num_steps-1;
    arma::mat A = tri_matrix(N, false);
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int iterations;
    int maxiter;
    bool converged;

    int k;
    int l;

    jacobi_eigensolver(A, pow(10,-20), eigenvalues, eigenvectors, maxiter, iterations, converged);

    //std::cout<< iterations; 
    //std::cout<< D;
    // //std::cout<< A;
    //std::cout<< "#################################################################################################################################################";
    //std::cout<< eigenvalues;
    //std::cout<< eigenvectors;

    //std::cout<< "#################################################################################################################################################";

    std::cout << eigenvalues;
    // arma::mat eigenvectrue;
    // arma::eig_sym(eigenvaltrue, eigenvectrue, A);
    // std::cout << eigenvectrue;
    // std::cout << eigenvaltrue;

    //std::cout<< "#################################################################################";
    //std::cout<< A;



        // Set a filename
    eigenvalues.save("Task6_output_eigenvalues100");

    eigenvectors.save("Task6_output_eigenvectors100");

    return 0;
}
