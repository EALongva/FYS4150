// test file project 5

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/rossby.hpp"

int main(){

  // testing rossby class
  rossby ross(3);
  ross.function();

  // generating an armadillo array to test data saving and plotting
  arma::vec r(10,arma::fill::randu);
  r.save("results/data/example.csv", arma::csv_ascii);

  // test the algorithms for calculating the laplacian from project 1

  return 0;
}
