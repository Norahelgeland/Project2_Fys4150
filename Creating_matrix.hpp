
#include <armadillo>
#include <iostream>
#include <math.h>



arma::mat tri_matrix(double N, bool dense);


arma::mat tri_matrix(double N, bool dense){

    arma::mat A = arma::mat(N,N).fill(0.);

    if(!dense){

        //Create tridiagonal matrix
        double h = 1/N;   //Our step_size
        double a = -1/pow(h,2);
        double b = 2/pow(h,2);
 
        A(0,0)=b;
        A(0,1)=a;
        A(N-1, N-1)=b;
        A(N-1, N-2)=a;

        for(int i=1; i < N-1; i++){
            for(int j=1;j<N-1; j++){

                if(i==j){
                    A(i,j)=b;
                    A(i,j-1)=a;
                    A(i,j+1)=a;
                
                }      
            }
        }
    }

    else{
        // Generate random N*N matrix
        //arma::mat 
        A = arma::mat(N, N).randn();  

        // Symmetrize the matrix by reflecting the upper triangle to lower triangle
        A = arma::symmatu(A);  

        }

    return A;
}



