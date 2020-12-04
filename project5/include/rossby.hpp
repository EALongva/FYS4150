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
  arma::vec psi_init;

  // streamfunction time evolution matrix (1D, x - spacial, y - time evol)
  arma::mat psi;

  // Initializer
  rossby(arma::vec psi_in);
  rossby(double u_in, double v_in);
  rossby(double u_in, double v_in, arma::vec psi_in);


  // Functions
  void function();
  void timestep_forward(double &zeta_forward, double psi_forward, double psi_backward, double deltat, double deltax)
  void timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward, double deltat, double deltax)
};

#endif
