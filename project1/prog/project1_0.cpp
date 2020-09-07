// project 1 algo

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
inline double exact(double x) {return 1.0-(1.0-std::exp(-10))*x-std::exp(-10*x);}

arma::vec tridiag_simp(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n, arma::fill::zeros);
  // updating the diagonal vector elements "forward sweep"
  c(0) = c(0)/b(0);
  d(0) = d(0)/b(0);
  for (int i = 1; i <= n-2; i++){
    c(i) = c(i)/(b(i)-a(i-1)*c(i-1));
    d(i) = (d(i)-a(i-1)*d(i-1))/(b(i)-a(i-1)*c(i-1));
  }
  //backward substitution computes the resulting vector v
  v(n-1) = (d(n-1)-a(n-2)*d(n-2))/(b(n-1)-a(n-2)*c(n-2)); //v(n-1) = d(n-1)
  for (int i = n-2; i >= 0; i--) v(i) = d(i)-c(i)*v(i+1);
  return v;
}

/*arma::vec tridiag_spec(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n);
  v.zeros();
  // updating the diagonal vector elements "forward sweep" (not looping over d(n-1))
  c(0) = c(0)/b(0);
  d(0) = d(0)/b(0);
  for (int i = 1; i <= n-2; i++){
    c(i) = c(i)/(b(i)-a(i-1)*c(i-1));
    d(i) = (d(i)-a(i-1)*d(i-1))/(b(i)-a(i-1)*c(i-1));
  }
  //backward substitution computes the resulting vector v

  v(n-1) = d(n-1); // this is incorrect
  for (int i = n-2; i >= 0; i--) v(i) = d(i)-c(i)*v(i+1);
  return v;
}*/

int main(int argc, char* argv[])
{
  // doing some clock work
  clock_t start;
  clock_t finish;
  start = clock();

  int n = 100; //length of arrays
  double h = 1/(n+1.0); //step length
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
  for (int i = 1; i < n+1; i++){
    d(i-1) = hh * f(i*h);
  }
  arma::vec v = tridiag_simp(a,b,c,d,n);
  //std::cout << v << std::endl;

  //ending the clock work
  finish = clock();
  std::cout << std::setprecision(32) << " timeused = " << ((double) (finish - start)/CLOCKS_PER_SEC) << std::endl;

  // writing to file
  std::string filename = argv[1];
  ofile.open(filename.c_str());
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << "number of steps," << "time," << std::endl;
  ofile << n << ", " << ((double) (finish - start)/CLOCKS_PER_SEC) << std::endl;
  ofile << "Numeric Solution," << "Exact Solution," << "Relative Error,"  << std::endl;
  for (int i = 1; i < n+1; i++){
    double RelError = std::fabs( (exact(i*h)-v(i-1))/exact(i*h));
    ofile << std::setprecision(12) << v(i-1) << ",";
    ofile << std::seatprecision(12) << exact(i*h) << ",";
    ofile << std::setprecision(12) << RelError << std::endl;
  }
  ofile.close();

  return 0;
}
