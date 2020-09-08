// project 1 tridiagonal special case algo, LINK -larmadillo to get LU working

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

arma::vec tridiag_LU(arma::mat A, arma::vec d, int n){

  // setting up the v vector, giving it the same length as input vector b
  // v is the solution vector (u in project description)
  arma::vec v(n, arma::fill::zeros);

  // also setup the intermidiate step vector y to solve the LU problem
  arma::vec y(n, arma::fill::zeros);

  // using the LU decomp from armadillo
  arma::mat L;
  arma::mat U;
  arma::mat P;
  arma::lu(L, U, P, A);

  // using arma::solve to solve the linear problem (2 steps)
  d = P * d; // unnecessary for our matrix ?
  y = arma::solve(L,d);
  v = arma::solve(U,y);
  
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

    // doing some clock work
    clock_t start;
    clock_t finish;
    start = clock();

    int n = ((int) std::pow(10.0, i)); //length of arrays
    double h = 1/(n+1.0); //step length
    double hh = h*h; // h^2

    /*
    // filling all the tremendously beautiful arrays with the best values
    arma::vec a(n-1, arma::fill::zeros);
    a = a - 1.0;
    arma::vec b(n, arma::fill::zeros);
    b = b + 2.0;
    arma::vec c(n-1, arma::fill::zeros);
    c = c - 1.0;
    */

    // declaring the known vector d (solution to Ax = d)
    arma::vec d(n, arma::fill::zeros);

    // filling the d array with function values times h^2 from x=0 to x=1
    for (int j = 1; j < n+1; j++){
      d(j-1) = hh * f(j*h);
    }
    // lets attempt to create our big and gorgeous matrix
    arma::mat A(n,n,arma::fill::zeros);
    // filling the diagonal elements
    A.diag(0) += 2.0; A.diag(-1) += -1.0; A.diag(1) += -1.0;
    //A.print();

    //computing the solution, given as vector v
    arma::vec v = tridiag_LU(A,d,n);

    //ending the clock work
    finish = clock();
    //std::cout << std::setprecision(32) << " timeused = " << ((double) (finish - start)/CLOCKS_PER_SEC) << std::endl;


    // writing to file

    std::string fileout = filename;
    std::string argument = std::to_string(i);
    fileout.append(argument);
    fileout.append(".csv");
    ofile.open(fileout.c_str()); //  .c_str()
    ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
    ofile << "steps," << "time," << std::endl;
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
