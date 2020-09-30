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

int main(){

  int N = 100;

  arma::mat A; // initial matrix A

  A = BBmatrix(N); // obtaining the initial buckling beam matrix

  arma::vec arma_eig_val;
  arma::mat arma_eig_vec;

  // keep in mind that eig_sym returns eigenvalues sorted
  arma::eig_sym(arma_eig_val, arma_eig_vec, A);

  arma::vec eigvals(N, arma::fill::zeros);
  arma::vec eigvec(N, arma::fill::zeros);

  BBanalyticEig(N, eigvals, eigvec);

  // writing to file
  arma::mat M(4,N,arma::fill::zeros);
  std::vector<std::string> v = {"analytic eigval", "arma eigval", "analytic eigvec", "arma eigvec"};
  std::string filename = "dat/analyticalEigComparisonBB.csv";

  M.row(0) = eigvals.t();
  M.row(1) = arma_eig_val.t();
  M.row(2) = eigvec.t();
  M.row(3) = arma_eig_vec.col(0).t();

  ToFile(M, v, filename);

  return 0;
}
