// functon library for project 2

#include <iostream>
#include <armadillo>
#include "../include/utils.hpp"

arma::mat maketridiag(double a, double b, double c, int N){

  arma::mat A(N,N,arma::fill::zeros);

  A.diag(1) += a;   // upper off diagonal
  A.diag(0) += b;   // diagonal
  A.diag(-1) += c;  // lower off diagonal

  return A;

}

void maximum_indices(arma::mat A, int N, int& k, int& l){

  double max1, max2, max; // declaring variables
                          // 1 refers to upper triangular matrix, 2 lower tri.mat.
  max1 = 0.0; max2 = 0.0;

  int k1, l1, k2, l2; // indices, 1 upper, 2 lower
  k1 = 0; k2 = 0; l1 = 0; l2 = 0; k = 0; l = 0;

  for(   int i = 0;   i < N-1;  i++){ // INDEXING -> A(i,j) = A(k,l)
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

arma::mat rotation(arma::mat A, int N, int k, int l){

  // Peforming rotation on maximum non-diagonal element

  double a_ik, a_il, a_kk, a_ll, a_kl;
  a_kk = A(k,k); a_ll = A(l,l); a_kl = A(k,l);

  double tau, t, c, s, cc, ss, cs;
  tau = (a_ll - a_kk)/(2*a_kl);

  if( tau > 0 ){
    t = 1./(tau + sqrt(1 + tau*tau));
  }

  else{
    t = -1./(-tau + sqrt(1 + tau*tau));
  }
  //t = tau + sqrt(1 + tau*tau);

  c = 1./sqrt(1 + t*t);
  s = t*c;
  cc = c*c; ss = s*s; cs = c*s;


  for(int i=0;i < N;i++){

    if(i != k && i != l){

      a_ik = A(i,k); a_il = A(i,l);
      A(i,k) = a_ik*c - a_il*s;
      A(i,l) = a_il*c + a_ik*s;
      A(k,i) = A(i,k);
      A(l,i) = A(i,l);

    }
  }

  A(k,k) = a_kk*cc - 2*a_kl*cs + a_ll*ss;
  A(l,l) = a_ll*cc + 2*a_kl*cs + a_kk*ss;
  A(k,l) = 0.0;
  A(l,k) = 0.0;

  return A;

}

arma::mat jacobimethod(arma::mat A, int N, int eps){

  double epsilon = std::pow(10.,eps); // taking the eps argument as a power
  int iterations = 0; // iteration counter
  double max = 10.0; // placeholder value
  int k, l; // declaring the indices k and l

  while (max > epsilon){

    iterations ++;

    // Finding the maximum non-diagonal element
    maximum_indices(A,N,k,l); // retrieves the maximum indices k,l
    max = std::fabs(A(k,l)); // calculating the new maximum off-diag element

    //std::cout << max << std::endl;

    // Peforming rotation on maximum non-diagonal element
    A = rotation(A, N, k, l);

    //std::cout << iterations << std::endl;

  }

  return A;

}
