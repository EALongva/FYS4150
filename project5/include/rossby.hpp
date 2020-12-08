#ifndef ROSSBY_H
#define ROSSBY_H

#include <iostream>
#include <cmath>
#include <armadillo>

class rossby
{
public:
  // Properties
  double deltax;
  double deltat;
  int xdim;
  int tdim;

  // // initial velocities
  // double u, v;

  // // initial psi function
  // arma::vec psi_init;

  // streamfunction time evolution matrix (1D, x - spacial, y - time evol)
  arma::mat psi;

  // vorticity time evolution matrix (1D, x - spacial, y - time evol)
  arma::mat zeta;

  // Initializer
  rossby(double dx, double dt, int J, int N);
  // rossby(double u_in, double v_in);
  // rossby(double u_in, double v_in, arma::vec psi_in);


  // Functions
  void function();
  void initialize_wave(bool sineWave, double sigma, double x0);
  void precalculate_offdiag(arma::vec &c_new);
  void gaussian_elimination(arma::vec &psi_forward, arma::vec zeta, arma::vec c_new);
  void timestep_forward(double &zeta_forward, double psi_forward, double psi_backward);
  void timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward);
};

#endif
