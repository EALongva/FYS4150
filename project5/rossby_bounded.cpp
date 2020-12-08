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
  rossby ross(0.5, 0.5, 10, 10);
  ross.initialize_wave(true, 0, 0);
  
  cout << ross.psi.col(0) << endl;

  return 0;
}
