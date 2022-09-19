#include <math.h>
#include <armadillo>
#include <iostream>
#include "Creating_matrix.hpp"
#include "Jacobi.hpp"



int main(){

    double num_size = 100;
    arma::vec N_vec = arma::vec(num_size-2);
    arma::vec num_iter = arma::vec(num_size-2);

    for(double i=2; i < num_size; i++){

        double N = i;

        arma::vec eigenvalues;
        arma::mat eigenvectors;
        int iterations;
        int maxiter;
        bool converged;
        int k;
        int l;
        arma::mat A = tri_matrix(N, true);
        jacobi_eigensolver(A, pow(10,-20), eigenvalues, eigenvectors, maxiter, iterations, converged);
        num_iter(N-2)=iterations;
        N_vec(N-2)=N;
        
        

    }


   
        // Set a filename
        std::string filename = "Task5_output_dense.txt";

        // Create and open the output file. Or, technically, create 
        // an "output file stream" (type std::ofstream) and connect 
        // it to our filename.
        std::ofstream ofile;
        ofile.open(filename);

        int width = 12;
        int prec  = 4;

        // Loop over steps
        for (int i = 0; i < N_vec.size(); i++)
        {
        // Write a line with the current x and y values (nicely formatted) to file
        ofile << std::setw(width) << std::setprecision(prec) << std::scientific << N_vec[i]
                << std::setw(width) << std::setprecision(prec) << std::scientific << num_iter[i]
                << std::endl;
        }  
        ofile.close();

    return 0;

}



