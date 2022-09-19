#include <armadillo>
#include <iostream>
#include <math.h>

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);


double max_offdiag_symmetric(const arma::mat& A, int& k, int& l){

    double max_value = 0;

    for(int i=0; i < A.n_rows-1; i++){
        for(int j=0;j < A.n_cols-1; j++){

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

