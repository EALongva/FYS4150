// example file project 4

#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
//#include "include/exampleClass.hpp"
//#include "include/utils.hpp"
#include "include/ising.hpp"

int main(){
  //srand(42);
  //srand(time(NULL));
  /*

  // testing utils
  example();

  // testing exampleClass
  exampleClass test(3);
  test.function();

  // generating an armadillo array to test data saving and plotting
  arma::vec r(10,arma::fill::randu);
  r.save("results/data/example.csv", arma::csv_ascii);
  */

  // generating random numbers using m
  /*std::mt19937 rng(42);
  std::uniform_real_distribution<double> p(0,1);
  for (int i=0; i<10; i++) std::cout << p(rng) << std::endl;
  */


  // testing metro class
  int N_in      = 2;
  int seed      = 42;
  double T_in   = 1.0;
  double J_in   = 1.0;

  Ising M(N_in, T_in, J_in, seed);


  M.genState();
  M.energy();

  std::cout << rand() << std::endl;

  return 0;
}









//
