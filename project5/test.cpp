// test file project 5

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/rossby.hpp"

int main(){

  //  testing rossby class
  //ross rossby();
  //ross.function();

  //  generating an armadillo array to test data saving and plotting
  //arma::vec r(10,arma::fill::randu);
  //r.save("results/data/example.csv", arma::csv_ascii);

  //  test the algorithms for calculating the laplacian from project 1
  double dx = 0.03; //spacial step length
  double x0 = 0;
  double x1 = 1;
  double xn_temp = ((x1 - x0) / dx);
  int xn = (int) xn_temp + 0.5;

  std::cout << "spacial steps:  " << xn << std::endl;

  int j;
  for (double i=x0; i<x1; i+=dx){
    j += 1;
  }

  std::cout << "forward loop steps:  " << j << std::endl;

  //arma::vec testpsi = arma::linspace()

  return 0;
}
