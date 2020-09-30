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

  // Defining the elements of the tridiagonal matrix
  vec d(N, fill::zeros); // Vector of diagonal elements
  vec a(N, fill::zeros); // Vector of off-diagonal elements
  vec aa(N, fill::zeros); // Vector of off-diagonal elements squared
  for (int i = 0; i < N; i++){
    d(i) = 2;
    a(i) = -1;
    aa(i) = a(i)*a(i);
  }
  a(0) = aa(0) = 0; // Setting the first elements of a and a^2 to zero
  double eps, nonzero;
  int m1,m2;
  eps = pow (10.,-10); // The precision to which we want to know our lambda
  nonzero = pow (10.,-308); // The smallest number on the computer not zero
  m1 = 0; // The index of the lowest eigenvalue we want to find
  m2 = N-1; // The index of the largest eigenvalue we want to find

  // Setting up vectors to store information
  vec x(N, fill:: zeros); // Vector to store eigenvalues/upper bounds
  vec lb(N, fill::zeros); // Vector to store lower bounds for eigenvalues

  // Setting initial search interval xmin and xmax
  double xmax, xmin;
  xmin = d(N-1) - abs(a(N-1));
  xmax = d(N-1) + abs(a(N-1));
  for (int i = N-2; i > 0; i--){
    double h = abs(a(i)) + abs(a(i+1));
    if (d(i) + h > xmax){
      xmax = d(i) + h;
    }
    if (d[i] - h < xmax){
      xmin = d(i) - h;
    }
  }
  // Setting lower bound and upper bound vectors equal to xmin and xmax
  for (int i = m1; i < m2+1; i++){
    x(i) = xmax;
    lb(i) = xmin;
  }

  double xl; // Declaring lower bound x, xl
  double xu; // Declaring upper bound x, xu
  double xm; // Declaring middle x, xm
  double r; // Declaring ratio between polynomials
  int z, c; // Declaring counter of iterations z and counter of roots c

  xu = xmax;
  z = 0;

  clock_t start, finish; // declaring start and finish time
  start = clock(); // start time of computing

  for (int k = m2; k > m1-1; k--){ // Looping over all eigenvalues k

    xl = xmin;

    for (int i = k; i > m1-1; i--){
      if (xl < lb(i)){
        xl = lb(i); // Setting the lower bound to the largest lower bound
      }
      if (xu > x(k)){
        xu = x(k); // Setting the upper bound to the smallest upper bound
      }
    }

    while (xu - xl > eps && z < (int) pow(10,5)){

      xm = (xl + xu)/2.; // Setting middle point
      z += 1;

      // Sturm sequence
      c = 0; r = 1; // We start with r = 1 since P_0 = 1
      for (int i = 0; i < N; i++){
        if (r != 0){
          r = d(i) - xm - aa(i)/r;
        }
        else{
          r = d(i) - xm - aa(i)/nonzero; // avoid division by zero
        }
        if (r < 0){
          c += 1; // Counting roots
        }
      }

      // If the k-th root is larger than xm, set lower bound to xm
      if (c < k + 1){
        if (c < m1 + 1){  // If there are less roots than lowest number of roots we want to find
          xl = lb(m1) = xm; // - make xm the lowest bound for all eigenvalues we want to find
        }
        else{ // If not
          xl = lb(c) = xm; // - make xm the lower bound for the eigenvalue with index larger than the number of roots
          if (x(c-1) > xm){
            x(c-1) = xm; // make xm the upper bound for the eigenvalue with index lower than the number of roots
          }
        }
      }

      // If the k-th root is smaller than xm, et upper bound to xm
      else{
        xu = xm;
      }

    }
    x(k) = (xl + xu)/2.; // Set k-th eigenvalue (or the (k-1)-th upper bound) to xm
  }
  finish = clock(); // finish time of computing
  double time = ((double)(finish - start)/CLOCKS_PER_SEC);
  cout << setprecision(32) << "Time used for " << N << " gridpoints: "<< time << endl;
  x.print();
  lb.print();
  cout << z << endl;
}
