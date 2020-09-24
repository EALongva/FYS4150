// functon library for project 2

#include <iostream>
#include <armadillo>
#include "../include/utils.hpp"

arma::mat maketridiag(double a, double b, double c, int N){

  arma::mat A(N,N,arma::fill::zeros);

  A.diag(-1) += a;  // lower off diagonal
  A.diag(0) += b;   // diagonal
  A.diag(1) += c;   // upper off diagonal

  return A;

}

void maximum_indices(arma::mat A, int N, int& k, int& l){

  double max1, max2, max; // declaring variables
                          // 1 refers to upper triangular matrix, 2 lower tri.mat.
  max1 = 0.0; max2 = 0.0;

  int k1, l1, k2, l2; // indices, 1 upper, 2 lower
  k1 = 0; k2 = 0; l1 = 0; l2 = 0; k = 0; l = 0;

  for(   int i = 0;   i < N-1;  i++){ //
    for( int j = i+1; j < N;    j++){ // Finding the maximum of the upper non-diagonal

      if ( fabs(A(i,j) ) > max1 ){

        max1 = fabs(A(i,j));
        std::cout << max1 << std::endl;
        k1 = i; l1 = j;
        std::cout << k1 << ", " << l1 << std::endl;
      }
    }
  }

  for(   int i = 1; i < N; i++){
    for( int j = 0; j < i; j++){ // Finding the maximum of the lower non-diagonal

      if ( fabs(A(i,j)) > max2 ){

        max2 = fabs(A(i,j));
        std::cout << max2 << std::endl;
        k2 = i; l2 = j;
        std::cout << k2 << ", " << l2 << std::endl;
      }
    }
  }

  if(max1 > max2){ // checking which maximum value (upper (1) or lower (2) tri.matrix)
                   // is greatest, max is then filled with the final maximum value
    max = max1;
    k = k1; l = l1;
  }

  else{
    max = max2;
    k = k2; l = l2;
  }
  std::cout << k << ", " << l << std::endl;
}
