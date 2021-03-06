// trying to figure out how to structure and solve exc 2b
//
// How to compile:
//          >>g++ -std=c++17 -c src/utils.cpp [creates utils.o in main dir]
//          >>g++ -std=c++17 -c main.cpp      [creates main.o in main dir]
//          >>g++ -std=c++17 utils.o main.o -o program.exe -l armadillo
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
#include "time.h"
#include "include/utils.hpp"

int main(int argc,char* argv[]){
  if (argc < 2){
    std::cout << "not enough arguments" << std::endl;
    exit(1);
  }

  int N = std::atoi(argv[1]);

  /*
  float rho0 = 0.0;
  float rhoN = 1.0;
  float h = (rhoN - rho0)/((double) N);
  float hh = h*h;
  //std::cout << h << hh << std::endl;

  float d = 2/hh;
  float a = -1/hh;
  */

  //arma::mat A = maketridiag(a, d, a, N);
  //A.print();

  //arma::mat D = arma::eig_sym(A);
  //D.print();

  /*
  arma::mat HO;
  HO = HOmatrix(rho0, rhoN, N);
  HO.print();
  */

  // testing the ToFileFunction

  /*
  arma::vec v1 = HO.diag(0);
  arma::vec v2 = HO.diag(1);
  arma::vec v3 = HO.diag(-1);
  */


  int eps = -8;
  int N_min = 5;
  int N_step = 10;
  int N_points = 13;
  std::string filename = "dat/SimTransCountN_13.csv";

  SimTransCount(eps, N_min, N_step, N_points, filename);



  /*
  arma::mat A(3, N, arma::fill::randu);

  std::string filename;
  filename = "data.csv";

  std::vector<std::string> v = {"n1", "n2", "n3"};

  ToFile(A, v, filename);

  arma::mat A = maketridiag(a,d,a,N);
  //A.print();


  //int k;
  //int l;
  //maximum_indices(A, N, k, l);
  //std::cout << k << ", " << l << std::endl;

  // timing!
  clock_t start;
  clock_t finish;

  int eps = -8;
  int iterations;
  start = clock();
  arma::mat B = jacobimethod(A,N,eps,iterations);
  finish = clock();
  //B.print();

  // timing end!
  double T = ((double) (finish - start)/CLOCKS_PER_SEC );

  std::cout << "iterations: " << iterations << "  time elapsed: " << T << std::endl;
  */


}
