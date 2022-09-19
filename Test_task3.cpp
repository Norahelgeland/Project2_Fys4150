
#include <armadillo>
#include <iostream>
#include <math.h>


int main(){

    arma::mat T = arma::mat(4,4).fill(0);
    T.diag() += 1;

    T(0,3) = 0.5;
    T(3,0) = 0.5;
    T(2,1) = -0.7;
    T(1,2) = -0.7;

    T.save("fileT");

    //int k;
    //int l;

    //double max_value = max_offdiag_symmetric(T, k, l);

}