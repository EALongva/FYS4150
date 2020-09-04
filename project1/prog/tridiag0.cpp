// lol

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include "time.h"

// object for output files
std::ofstream ofile;

inline double f(double x){return 100.0*std::exp(-10*x);}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

arma::vec tridiag_simp(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n);
  v.zeros();
  // updating the diagonal vector elements "forward sweep"

  for (int i = 1; i <= n-2; i++){
    double w = a(i) / b(i-1);
    b(i) = b(i) - w*c(i-1);
    d(i) = d(i) - w*d(i-1);
  }
  //backward substitution computes the resulting vector v

  v(n-1) = d(n-1) / b(n-1);
  for (int i = n-2; i >= 0; i--) v(i) = (d(i)-c(i)*v(i+1))/b(i);
  return v;
}

int main(int argc, char* argv[])
{
  // doing some clock work
  clock_t start;
  clock_t finish;
  start = clock();

  int n = 100; //length of arrays
  double h = 1/(n-1.0); //step length squared
  double hh = h*h; // h^2
  //std::cout << hh << std::endl;
  // filling all the tremendously beautiful arrays with the best values
  arma::vec a(n-1, arma::fill::zeros);
  a = a - 1.0;
  arma::vec b(n, arma::fill::zeros);
  b = b + 2.0;
  arma::vec c(n-1, arma::fill::zeros);
  c = c - 1.0;
  arma::vec d(n, arma::fill::zeros);
  // filling the d array with function values times h^2 from x=0 to x=1
  for (int i = 0; i < n; i++){
    d(i) = hh * f(i*h);
  }
  arma::vec v = tridiag_simp(a,b,c,d,n);
  //std::cout << v << std::endl;

  //ending the clock work
  finish = clock();
  std::cout << " timeused = " << ((double) (finish - start)/CLOCKS_PER_SEC) << std::endl;

  // writing to file
  std::string filename = argv[1];
  ofile.open(filename.c_str());
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << "numSolution    " << "exactSolution    " << "relativeError    "  << std::endl;
  for (int i = 0; i < n; i++){
    double RelError = std::fabs( (exact(i*h)-v(i))/exact(i*h) );
    ofile << std::setw(15) << std::setprecision(10) << v(i);
    ofile << std::setw(15) << std::setprecision(10) << exact(i*h);
    ofile << std::setw(15) << std::setprecision(10) << RelError << std::endl;
  }
  ofile.close();

  return 0;
}
