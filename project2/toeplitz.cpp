#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>
#include "time.h"

// use namespace for output and input
using namespace std;
using namespace arma;
// object for output files
ofstream ofile;

int main(int argc, char *argv[]){
  int N;
  string ofilename;

  if(argc <= 1){
    cout << "Two few arguments in " << argv[0] <<
    ". Read in name of output file and the number of mesh points." << endl;
    exit(1);
  }
  else{
    ofilename = argv[1]; // Output file name
    N = atoi(argv[2]); // Matrix dimension
  }

  // Defining fixed values
  double rho0, rhoN, h, hh, d, a;
  rho0 = 0.0;         // Minimum rho value
  rhoN = 1.0;         // Maximum rho value
  h = (rho0-rhoN)/N;  // Step size
  hh = h*h;
  d = 2./hh;        // Diagonal element constant
  a = -1./hh;       // Non-diagonal element constant

  // Setting up tridiagonal Toeplitz matrix using armadillo
  mat Toeplitz = zeros<mat>(N,N);
  Toeplitz(0,0) = d;
  Toeplitz(0,1) = a;
  for(int i = 1; i < N-1;i++){
    Toeplitz(i,i-1) = a;
    Toeplitz(i,i) = d;
    Toeplitz(i,i+1) = a;
  }
  Toeplitz(N-1,N-2) = a;
  Toeplitz(N-1,N-1) = d;

  // Diagonalizing the matrix with armadillo
  vec eigenvalues(N);
  eig_sym(eigenvalues,Toeplitz);
  double pi = 2*acos(0.0);
  // Opening file and writing out results
  ofile.open(ofilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << "Results from Diagonalization of Toeplitz Matrix" << endl;
  ofile << "         exact:     armadillo:         error:" << endl;
  for (int i = 1; i < N+1; i++){
    double exact = d + 2*a*cos(i*pi/(N+1));
    double Error = fabs(eigenvalues(i-1)-exact); // Compute relative error
    ofile << setw(15) << setprecision(8) << exact;
    ofile << setw(15) << setprecision(8) << eigenvalues(i-1);
    ofile << setw(15) << setprecision(8) << Error << endl;
  }
  ofile.close();

  // Jacobi's method
  double epsilon = pow (10.,-8);

  // Finding the maximum non-diagonal element
  double max1, max2, max;
  int k1, l1, k2, l2, k, l;
  max1 = 0.0; max2 = 0.0;
  for(int i = 0;i < N;i++){
    for(int j = i+1;j < N;j++){ // Finding the maximum of the upper non-diagonal
      if (fabs(Toeplitz(i,j))>max1){
        max1 = fabs(Toeplitz(i,j));
        k1 = i; l1 = j;
  }}}
  for(int i = 1;i < N;j++){
    for(int j = 0;j < i;j++){ // Finding the maximum of the lower non-diagonal
      if fabs(Toeplitz(i,j))>max2){
        max2 = fabs(Toeplitz(i,j));
        k2 = i; l2 = j;
  }}}
  if(max1 > max2){
    max = max1;
    k = k1; l = l1;
  }
  else{
    max = max2;
    k = k2; l = l2;
  }

// Rotation
double a_kl, a_kk, a_ll, tau;
a_kl = Toeplitz(k,l);
a_kk = Toeplitz(k,k);
a_ll = Toeplitz(l,l):
tau = (a_ll - a_kk)/(2*a_kl);


}
