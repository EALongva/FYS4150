// example file project 4

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/exampleClass.hpp"
#include "include/utils.hpp"

int main(){

  // testing utils
  example();

  // testing exampleClass
  exampleClass test(3);
  test.function();

  // generating an armadillo array to test data saving and plotting
  arma::vec r(10,arma::fill::randu);
  r.save("results/data/example.csv", arma::csv_ascii);

  return 0;
}
