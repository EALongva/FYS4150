/* script solving project 2 b) Finding the similarity transformations needed
to diagonalize the buckling beam matrix */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <cmath>
#include <vector>
#include "time.h"
#include "../include/utils.hpp"

int main(){

  int eps = -8;
  int N_min = 5;
  int N_step = 10;
  int N_points = 13;
  std::string filename = "dat/SimTransCount.csv";

  SimTransCount(eps, N_min, N_step, N_points, filename);

  return 0;
}
