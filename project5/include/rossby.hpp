#ifndef ROSSBY_H
#define ROSSBY_H

#include <iostream>
#include <cmath>
#include <armadillo>

class rossby
{
public:
  // Properties
  int number;

  // initial velocities
  double u, v;

  // initial psi function
  arma::vec psi;

  // Initializer
  rossby(arma::vec psi_in);
  rossby(double u_in, double v_in);
  rossby(double u_in, double v_in, arma::vec psi_in);
  

  // Functions
  void function();
};

#endif
