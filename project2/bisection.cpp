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
  mat A = zeros<mat>(N,N);
  A(0,0) = d;
  A(0,1) = a;
  for(int i = 1; i < N-1;i++){
    A(i,i-1) = a;
    A(i,i) = d;
    A(i,i+1) = a;
  }
  A(N-1,N-2) = a;
  A(N-1,N-1) = d;

  
}
