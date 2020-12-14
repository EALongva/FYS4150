#ifndef ROSSBY_H
#define ROSSBY_H

#include <iostream>
#include <cmath>
#include <armadillo>

class rossby
{
public:
  // Properties
  double endpos; // Length of spatial dimension
  double endtime; // Length of time period
  double deltax; // Grid space
  double deltat; // Time step
  int xdim; // Spatial dimension
  int tdim; // Temporal dimension
  arma::mat Psi; // Streamfunction matrix
  arma::mat Zeta; // Vorticity matrix

  // Initializer
  rossby(double dx, double dt, double endtime);

  // Functions
  void initialize_wave(bool sineWave, double sigma, double x0);
  void zeta_timestep_forward(double &zeta_forward, double zeta, double psi_forward, double psi_backward);
  void zeta_timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward);
  arma::vec precalculate_offdiag();
  void gaussian_elimination(int n, arma::vec c_new);
  void jacobis_method(int n, arma::vec zeta);
  void evolve_bounded(bool forwardStep);
  void evolve_periodic(bool forwardStep);
};

#endif
