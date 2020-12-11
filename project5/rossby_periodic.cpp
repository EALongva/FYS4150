#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/rossby.hpp"

using namespace arma;
using namespace std;

int main(){

  // testing rossby class
  rossby ross(0.02, 0.1, 10);
  ross.initialize_wave(true, 0, 0);
  ross.evolve_periodic(true);
  cout << ross.Zeta << endl;

  return 0;
}
