// functon library for project 2

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "time.h"
#include "../include/utils.hpp"

arma::mat maketridiag(double a, double b, double c, int N){

  arma::mat A(N,N,arma::fill::zeros);

  A.diag(1) += a;   // upper off diagonal
  A.diag(0) += b;   // diagonal
  A.diag(-1) += c;  // lower off diagonal

  return A;

}

arma::mat BBmatrix(int N){

  double rho0 = 0.0;    // static rho0 and rhoN
  double rhoN = 1.0;
  double h    = (rhoN - rho0)/((double) N);
  double hh   = h*h;

  double d    = 2/hh;
  double a    = -1/hh;

  arma::mat A(N,N,arma::fill::zeros);

  A.diag(0)   += d;   // diagonal
  A.diag(1)   += a;   // upper off diagonal
  A.diag(-1)  += a;   // lower off diagonal

  return A;

}

arma::mat HOmatrix(double rho0, double rhoN, int N){

  arma::mat HO(N, N, arma::fill::zeros);
  arma::vec rho = arma::linspace(rho0, rhoN, N);
  rho = arma::pow(rho,2); // gives us the rho^2_i values in an array 'rho'

  double h = (rhoN - rho0)/((double) N);
  double hh = h*h;

  double e = - 1/hh;
  double d = 2/hh;

  // filling in the values on the diagonals

  HO.diag(-1) += e;
  HO.diag(1) += e;

  HO.diag(0) = rho;
  HO.diag(0) += d;

  return HO;

}


void BBanalyticEig(int N, arma::vec& eigvals, arma::vec& eigvec){
  // returns the eigenvector for the lowest eigenvalue, rho (0,1)

  double h = 1/((double) N);
  double hh = h*h;
  double d = 2/hh;
  double a = -1/hh;
  double pi = M_PI;
  double arg = pi/((double) N);

  for (int j=1; j < N+1; j++){
    eigvals(j-1) = d + 2*a*std::cos(j*arg);
  }

  for (int j = 1; j < N+1; j++){
    eigvec(j-1) = std::sin(j*arg);
  }

  eigvec = eigvec / arma::norm(eigvec);




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
        //std::cout << max1 << std::endl;
        k = i; l = j;
        //std::cout << k1 << ", " << l1 << std::endl;
      }
    }
  }
/*
  for(   int i = 1; i < N; i++){
    for( int j = 0; j < i; j++){ // Finding the maximum of the lower non-diagonal

      if ( fabs(A(i,j)) > max2 ){

        max2 = fabs(A(i,j));
        //std::cout << max2 << std::endl;
        k2 = i; l2 = j;
        //std::cout << k2 << ", " << l2 << std::endl;
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
  */
  //std::cout << k << ", " << l << std::endl;
}

void rotation(arma::mat& A, arma::mat& R, int N, int k, int l){

  // Peforming rotation on maximum non-diagonal element
  //std::cout << "here" << std::endl;

  double a_ik, a_il, a_kk, a_ll, a_kl;
  a_kk = A(k,k); a_ll = A(l,l); a_kl = A(k,l);

  //std::cout << "here too" << std::endl;
  // elemnts of the eigenvector matrix R
  double r_ik, r_il;

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

  //std::cout << "here2" << std::endl;

  for(int i=0;i < N;i++){

    if(i != k && i != l){

      a_ik = A(i,k); a_il = A(i,l);
      A(i,k) = a_ik*c - a_il*s;
      A(i,l) = a_il*c + a_ik*s;
      A(k,i) = A(i,k);
      A(l,i) = A(i,l);

    }
    // eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);
    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;

  }

  A(k,k) = a_kk*cc - 2*a_kl*cs + a_ll*ss;
  A(l,l) = a_ll*cc + 2*a_kl*cs + a_kk*ss;
  A(k,l) = 0.0;
  A(l,k) = 0.0;

  //std::cout << "rotation fin" << std::endl;

}

arma::mat jacobimethod(arma::mat A, arma::mat& R, int N, int eps, int& iterations){

  double epsilon = std::pow(10.,eps); // taking the eps argument as a power
  // int iterations = 0; // iteration counter
  double max = 10.0; // placeholder value
  int k, l; // declaring the indices k and l

  arma::mat I(N, N, arma::fill::eye); // Identity matrix -> transforms to the
                                      // eigenvectors

  R = I;

  while (max > epsilon){

    iterations ++;

    // Finding the maximum non-diagonal element
    maximum_indices(A,N,k,l); // retrieves the maximum indices k,l
    max = std::fabs(A(k,l)); // calculating the new maximum off-diag element

    // Peforming rotation on maximum non-diagonal element
    rotation(A, R, N, k, l);

  }

  return A;

}

void ToFile(arma::mat A, std::vector<std::string> v, std::string filename){

  // Matrix A should have a shape of (n, N) with n as row elements, N as column
  // elements. The vector v should contain the labels for the values on each of
  // the rows in A.
  ///*
  if (A.n_rows != v.size()){
    std::cout << "shape of vector v does not match shape of matrix A" << std::endl;
    std::cout << "make sure elements in v match rows in matrix A" << std::endl;
    exit(1);
  }

  int n = v.size();
  int N = A.n_cols;

  // object for output files
  std::ofstream ofile;

  ofile.open(filename); //.c_str()

  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);

  for (int i = 0; i < n; i++){

    if(i == n-1){
      ofile << v[i]; // so that lhe rows do not end with a comma
    }

    else{
    ofile << v[i] << " ,  " ; // std::vector indexing v[] not arma convention
    }

  }

  ofile << std::endl;

  for (int j = 0; j < N; j++){ // this is a bit confusing, sorry
    for (int i=0; i < n; i++){

      if(i == n-1){
        ofile << std::setprecision(12) << A(i,j);
      }

      else{
      ofile << std::setprecision(12) << A(i,j) << ",";
      }

    }

    ofile << std::endl;

  }

  ofile.close();
  //*/

}

void SimTransCount(int eps, int N_min, int N_step, int N_points, std::string filename){
  // N_min: start N value, N_step: value added to N for each step,
  // N_points: number of times N_step is added to N_min (length
  // of resulting data file)

  arma::mat M(3, N_points, arma::fill::zeros); // matrix containing data to file
  arma::mat A, B, R;


  for (int i = 0; i < N_points; i++){

    int N = N_min + i*N_step;
    A = BBmatrix(N);

    clock_t start, finish;

    int iterations = 0;

    start = clock();
    B = jacobimethod(A, R, N, eps, iterations);
    finish = clock();

    double T = ((double) (finish - start)/CLOCKS_PER_SEC );

    // adding values for: N, iterations, T for this N to the M matrix
    M(0,i) = N;
    M(1,i) = iterations;
    M(2,i) = T;

  }

  std::vector<std::string> v = {"N", "iterations", "CPUTime"};
  ToFile(M, v, filename);

}

//
