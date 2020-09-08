// project 1 cpu time comparison of general and specialized (tridiag) algorithms

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
  arma::vec t(maxpower,arma::fill::zeros); //time array LU decomp solution

  for (int i = 1; i < maxpower+1; i++){

    int n = ((int) std::pow(10.0, i)); //length of arrays
    double h = 1/(n+1.0); //step length
    double hh = h*h; // h^2

    // known array A*v = d
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
    clock_t start;
    clock_t finish;

    start = clock();
    arma::vec v = tridiag_LU(A,d,n);
    finish = clock();

    // calculating the elapsed time during algorithm
    t(i-1) = ((double) (finish - start)/CLOCKS_PER_SEC);

  }
  ofile.open(filename.c_str());
  ofile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
  ofile << "steps," << "LU algo time," << std::endl;
  for (int i = 0; i < maxpower; i++){
    ofile << std::setprecision(12) << i+1 << ",";
    ofile << std::setprecision(12) << t(i) << std::endl;
  }
  ofile.close();
  return 0;
}
