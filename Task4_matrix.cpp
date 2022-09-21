

#include <armadillo>
#include <iostream>
#include <math.h>


int main(){

    arma::mat T = arma::mat(6,6).fill(0.);
    T.diag() += 1;

    T(0,5) = 0.5;
    T(5,0) = 0.5;

    T(4,1) = -0.7;
    T(3,2) = -0.7;
   
    T(2,3) = -0.7;
    T(1,4) = -0.7;
  

    T.save("file4");

}