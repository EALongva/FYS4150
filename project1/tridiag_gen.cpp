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

arma::vec tridiag_general(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n, arma::fill::zeros);
  // updating the diagonal vector elements "forward sweep"
  c(0) = c(0)/b(0);
  d(0) = d(0)/b(0); // 2 flops
  for (int i = 1; i <= n-2; i++){
    c(i) = c(i)/(b(i)-a(i-1)*c(i-1)); // 3 flops
    d(i) = (d(i)-a(i-1)*d(i-1))/(b(i)-a(i-1)*c(i-1)); // 5 flops
  }
  //backward substitution computes the resulting vector v
  v(n-1) = (d(n-1)-a(n-2)*d(n-2))/(b(n-1)-a(n-2)*c(n-2)); // 5 flops
  for (int i = n-2; i >= 0; i--) v(i) = d(i)-c(i)*v(i+1); // 2 flops
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

  for (int i = 1; i < maxpower+1; i++){

    int n = ((int) std::pow(10.0, i)); //length of arrays
    double h = 1/(n+1.0); //step length
    double hh = h*h; // h^2

    // filling all the tremendously beautiful arrays with the best values
    arma::vec a(n-1, arma::fill::zeros);
    a = a - 1.0;
    arma::vec b(n, arma::fill::zeros);
    b = b + 2.0;
    arma::vec c(n-1, arma::fill::zeros);
    c = c - 1.0;

    // known array A*v = d
    arma::vec d(n, arma::fill::zeros);
    // filling the d array with function values times h^2 from x=0 to x=1
    for (int j = 1; j < n+1; j++){
      d(j-1) = hh * f(j*h);
    }

    // initiating clock
    clock_t start;
    clock_t finish;
    start = clock();

    arma::vec v = tridiag_general(a,b,c,d,n);

    //ending the clock work
    finish = clock();
    std::cout << std::setprecision(32) << " timeused = " << ((double) (finish - start)/CLOCKS_PER_SEC) << std::endl;

    // writing to file
    std::string fileout = filename;
    std::string argument = std::to_string(i); // converting i to string
    fileout.append(argument); // appending the power to the filename
    fileout.append(".csv"); // adding the comma separated values extension
    ofile.open(fileout.c_str()); //  .c_str()
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "number of steps," << "time," << std::endl;
    ofile << n << ", " << ((double) (finish - start)/CLOCKS_PER_SEC) << std::endl;
    ofile << "Numeric Solution," << "Exact Solution," << "Relative Error,"  << std::endl;
    for (int i = 1; i < n+1; i++){
      double RelError = std::fabs( (exact(i*h)-v(i-1))/exact(i*h));
      ofile << std::setprecision(12) << v(i-1) << ",";
      ofile << std::setprecision(12) << exact(i*h) << ",";
      ofile << std::setprecision(12) << RelError << std::endl;
    }
    ofile.close();
  }
  return 0;
}
