/* script solving project 2 b) Finding the similarity transformations needed
to diagonalize the buckling beam matrix */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <vector>
#include "time.h"
#include "../include/utils.hpp"

int main(int argc,char* argv[]){

  if (argc < 3){
    std::cout << "Arguments required: N, rhoN" << std::endl;
    exit(1);

  }

  int       N = std::atoi(argv[1]);
  double rhoN = std::atof(argv[2]);
  int     eps = -8;

  // not sure how to get the eigenvectors but lets first find the eigenvalues

  arma::mat A; // initial matrix A
  arma::mat B; // diagonalized matrix B
  arma::mat R; //
  arma::vec u; // used to store the UNSORTED eigenvalues

  double rho0;
  rho0 = 0;
  A = HOmatrix(rho0, rhoN, N); // obtaining the initial qm harmonic oscillator matrix
  //A.print();

  int iterations; // not necessary for this problem
  B = jacobimethod(A, R, N, eps, iterations);

  u = B.diag();
  arma::vec eig_val = arma::sort(u,"ascend");

  // finding the index of the minimum eigenvalue (arma::index_min did not work)
  int i_min = 0;
  for (int i = 0; i < N-1; i++){
    if (u(i+1) < u(i_min)){
      i_min = i+1;
    }
  }
  //std::cout << i_min << " , " << u(i_min) << std::endl;

  arma::vec arma_eig_val;
  arma::mat arma_eig_vec;

  // keep in mind that eig_sym returns eigenvalues sorted
  arma::eig_sym(arma_eig_val, arma_eig_vec, A);

  // writing to file
  arma::mat M(4,N,arma::fill::zeros);
  std::vector<std::string> v = {"jacobi eigval", "arma eigval", "jacobi eigvec", "arma eigvec"};
  std::string filename = "dat/qm_eigen.csv";

  M.row(0) = eig_val.t();
  M.row(1) = arma_eig_val.t();
  M.row(2) = R.col(i_min).t();
  M.row(3) = arma_eig_vec.col(0).t();

  ToFile(M, v, filename);

  return 0;
}
