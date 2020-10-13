// class header for project 3

#ifndef SS_H
#define SS_H

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "time.h"
#include "../include/planet.hpp"

class ss
{
// calss variables
  public:
    int total_planets;
    double G;
    std::vector<planet> all_planets;
    std::vector<arma::mat> Xevo;
    std::vector<arma::mat> Vevo;

// constructor
  ss();

// functions
  void add(planet newplanet);
  void planet_names();
  void euler(double T, int N);
  void eulerStep(int I, int J, int i, double dt);
  void VVerlet(double T, int N);
  void test(int nr);

};




#endif
