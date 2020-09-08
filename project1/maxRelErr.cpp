// project 1 general tridiagonal algo

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

arma::vec tridiag_special(arma::vec b, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n, arma::fill::zeros);
  // updating the diagonal vector elements "forward sweep" (not looping over d(n-1))
  for (int i = 1; i < n+1; i++){
    b(i-1) = (i+1.0)/((double) i);
  }
  for (int i = 1; i < n; i++){
    d(i) = d(i) + d(i-1)/b(i-1);
  }
  //backward substitution computes the resulting vector v
  v(n-1) = d(n-1)/b(n-1);
  for (int i = n-2; i >= 0; i--) v(i) = (d(i) + v(i+1))/b(i);
  return v;
}

int main(int argc, char* argv[])
{
  if (argc <= 2){
    std::cout << "not enough arguments" << std::endl;
    std::cout << "arguments given: " << argc << "give filename and max power" << std::endl;
    exit(1);
  }

  std::string filename = argv[1];
  int maxpower = atoi(argv[2]);

  arma::vec relativeError(maxpower, arma::fill::zeros);

  for (int i = 1; i < maxpower+1; i++){

    int n = ((int) std::pow(10.0, i)); //length of arrays
    double h = 1/(n+1.0); //step length
    double hh = h*h; // h^2

    // filling all the tremendously beautiful arrays with the best values
    arma::vec b(n+1, arma::fill::zeros);
    b = b + 2.0;
    /*arma::vec a(n-1, arma::fill::zeros);
    a = a - 1.0;
    arma::vec c(n-1, arma::fill::zeros);
    c = c - 1.0;*/

    // known array A*v = d, and the closed form solution u
    arma::vec d(n, arma::fill::zeros);
    arma::vec u(n, arma::fill::zeros);

    // filling the d array with function values times h^2 from x=0 to x=1
    for (int j = 1; j < n+1; j++){
      d(j-1) = hh * f(j*h);
      u(j-1) = exact(j*h);

    }

    arma::vec v = tridiag_special(b,d,n); // computing special solution
    //arma::vec v = tridiag_general(a,b,c,d,n);

    // calculating the max relative error for each n
    arma::vec r = arma::abs((v-u)/u);
    double rmax = arma::max(r);
    relativeError(i-1) = std::log10(rmax);

  }

  ofile.open(filename.c_str());
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << "Steps," << "Relative Error" << std::endl;
  for (int i = 0; i < maxpower; i++){
    ofile << std::setprecision(12) << i+1 << ",";
    ofile << std::setprecision(12) << relativeError(i) << std::endl;
  }
  ofile.close();
  return 0;
}
