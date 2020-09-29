// trying to figure out how to structure and solve exc 2b
//
// How to compile:
//          >>g++ -std=c++17 -c src/utils.cpp [creates utils.o in main dir]
//          >>g++ -std=c++17 -c test.cpp      [creates test.o in main dir]
//          >>g++ -std=c++17 utils.o test.o -o program.exe -l armadillo
//
// note: -c flag creates the .o-files to link, -o flag let you name the program,
// -std choose c++ version, -l links the armadillo library

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <vector>
#include <cassert>
#include "time.h"
#include "include/utils.hpp"

void ASSERT_MAX_IND(){

  std::cout << "TESTING MAXIMUM_INDICES FUNCTION" << std::endl;
  std::cout << std::endl;

  int N = 5;
  double max = 10;
  int k, l;

  // testing positive max value
  double a, d;
  a = std::fabs(max);
  d = a + 1; // so that d > max

  arma::mat A;
  A = { {d,0,0,0,a},
        {0,d,0,0,0},
        {0,0,d,0,0},
        {0,0,0,d,0},
        {0,0,0,0,d}};

  std::cout << "matrix being tested: " << std::endl;
  A.print();


  maximum_indices(A,N,k,l);

  assert(A(k,l)==a);
  std::cout << "success, d = " << d << ", A(k,l) = " << a << std::endl;

  // testing positive max value

  a = std::fabs(max) * -1;
  d = a - 1; // so that d > max

  A = { {d,0,0,0,0},
        {0,d,0,0,0},
        {0,0,d,0,0},
        {0,0,0,d,0},
        {0,0,0,a,d}};

  std::cout << "matrix being tested: " << std::endl;
  A.print();


  maximum_indices(A,N,k,l);

  assert(A(k,l)==a);
  std::cout << "success, d = " << d << ", A(k,l) = " << a << std::endl;
  std::cout << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << std::endl;
}


void ASSERT_JACOBIMETHOD(){

  // test if our jacobimethod yields correct eigenvalues by comparison with
  // the eigenvalues obtained via armadillo's eig_sym function

  std::cout << "TESTING JACOBIMETHOD EIGENVALUES" << std::endl;
  std::cout << std::endl;

  int N = 5;
  int eps_power = -8;
  double eps = std::pow(10.,eps_power);
  int iterations;
  arma::mat A,B,R;

  B = { {2.0834,   0.8972,   0.9807,   1.8070,   1.0763},
        {0.8972,   0.6175,   0.7182,   0.9555,   0.3684},
        {0.9807,   0.7182,   0.9777,   1.1090,   0.5905},
        {1.8070,   0.9555,   1.1090,   1.7763,   0.9307},
        {1.0763,   0.3684,   0.5905,   0.9307,   1.0762}};
  std::cout << "matrix being tested:" << std::endl;
  B.print();
  std::cout << std::endl;

  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, B);

  std::cout << "armadillo eigenvalues: ";
  eigval.t().print();
  std::cout << std::endl;

  B = jacobimethod(B,R,N,eps_power,iterations);
  arma::vec jac_eigval;
  jac_eigval = B.diag();
  jac_eigval = arma::sort(jac_eigval);

  std::cout << "jacobi method eigenvalues: ";
  jac_eigval.t().print();
  std::cout << std::endl;

  arma::vec testvec;
  testvec = eigval - jac_eigval;
  double testvec_norm = arma::norm(testvec);

  assert(testvec_norm<eps);
  std::cout << "success, difference in eigen values: " << testvec_norm << std::endl;
  std::cout << "threshold used in jacobimethod: " << eps << std::endl;
  std::cout << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << std::endl;

}


int main(){

  ASSERT_MAX_IND();

  ASSERT_JACOBIMETHOD();

  std::cout << std::endl;
  std::cout << "All tasks failed successfully." << std::endl;

  return 0;
}
