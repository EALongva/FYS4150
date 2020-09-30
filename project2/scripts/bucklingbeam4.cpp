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

  int eps = -8;
  int Nmin = 10;
  int Nstep = 16;
  int Nmax = 160;
  int N;

  // not sure how to get the eigenvectors but lets first find the eigenvalues

  arma::mat A; // initial matrix A
  arma::mat B; // diagonalized matrix B
  arma::mat R; //
  arma::vec u; // used to store the UNSORTED eigenvalues
  arma::vec arma_eig_val;
  arma::mat arma_eig_vec;

  int iterations; // not necessary for this problem


  // CPU TIME COMPARISON jacobi vs armadillo eigsym
  clock_t staJac, finJac, staArm, finArm;

  arma::vec JacT(Nstep, arma::fill::zeros);
  arma::vec ArmaT(Nstep, arma::fill::zeros);

  for (int i=1; i<Nstep+1; i++){

    N = Nmin*i;
    A = BBmatrix(N); // obtaining the initial buckling beam matrix

    staJac = clock();
    B = jacobimethod(A, R, N, eps, iterations);
    finJac = clock();

    staArm = clock();
    arma::eig_sym(arma_eig_val, arma_eig_vec, A);
    finArm = clock();

    JacT(i-1) = ((double) (finJac - staJac)/CLOCKS_PER_SEC );
    ArmaT(i-1)= ((double) (finArm - staArm)/CLOCKS_PER_SEC );

  }

  u = arma::linspace(Nmin,Nmax,Nstep);

  // writing to file
  arma::mat M(3,Nstep,arma::fill::zeros);
  std::vector<std::string> v = {"N", "jacobi time", "arma time"};
  std::string filename = "dat/time.csv";

  M.row(0) = u.t();
  M.row(1) = JacT.t();
  M.row(2) = ArmaT.t();

  ToFile(M, v, filename);

  return 0;
}
