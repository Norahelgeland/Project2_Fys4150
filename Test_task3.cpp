#include <armadillo>
#include <iostream>
#include <math.h>
#include "Max_offdiag_symmetric.hpp"
#include <cassert>

void Test_function(); //Define the function

int main(){

    // Run the test code
    Test_function(); 
}

void Test_function(){

    // Define and create a T matrix 
    arma::mat T = arma::mat(4,4).fill(0);
    T.diag() += 1;

    T(0,3) = 0.5;
    T(3,0) = 0.5;
    T(2,1) = -0.7;
    T(1,2) = -0.7;

    T.save("TestMat"); // Saves matrix

    // Define values
    int k;
    int l;
    // Find the values trough the function
    double max_value = max_offdiag_symmetric(T, k, l);
    // Define the true values
    double true_max = 0.7;
    int true_k = 2;
    int true_l = 1;

    // Assertions to test if the found values are equal the true values
    assert(max_value=true_max);
    assert(k=true_k);
    assert(l=true_l);

}
