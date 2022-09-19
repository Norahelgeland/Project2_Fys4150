#include <armadillo>
#include <iostream>
#include <math.h>
#include "Creating_matrix.hpp"

int main(){
 
    double N = 6;

    arma::mat A = tri_matrix(N, false);

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