// trying to figure out how to structure and solve exc 2b
// run string:

#include <iostream>
#include <armadillo>
#include <cmath>
#include "include/utils.hpp"

int main(int argc,char* argv[]){
  if (argc < 2){
    std::cout << "not enough arguments" << std::endl;
    exit(1);
  }

  int N = std::atoi(argv[1]);
  float rho0 = 0.0;
  float rhoN = 1.0;
  float h = (rhoN - rho0)/((double) N);
  float hh = h*h;
  //std::cout << h << hh << std::endl;

  float d = 2/hh;
  float a = -1/hh;

  //arma::mat A = maketridiag(a, d, a, N);
  //A.print();

  //arma::mat D = arma::eig_sym(A);
  //D.print();

  arma::mat A(N,N,arma::fill::randu);
  A.print();

  int k;
  int l;
  maximum_indices(A, N, k, l);
  std::cout << k << ", " << l << std::endl;

}
