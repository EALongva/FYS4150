#ifndef ROSSBY_2D_H
#define ROSSBY_2D_H

#include <iostream>
#include <cmath>
#include <armadillo>

class rossby
{
public:
  // Properties
  double endposx; // Length of spatial domain in x-direction
  double endposy; // Length of spatial domain in y-direction
  double endtime; // Length of time period
  double deltax;  // Spatial x dimension
  double deltay;  // Spatial y dimension
  double deltat; // Temporal dimension
  int xdim; // x grid space
  int ydim; // y grid space
  int tdim; // time step
  arma::cube Psi; // Streamfunction matrix
  arma::cube Zeta; // Vorticity matrix

  // Initializer
  rossby(double dpos, double dt, double endtime);

  // Functions
  void initialize_wave(bool sineWave, double sigma, double x0, double y0);
  void zeta_timestep_forward(double &zeta_forward, double zeta, double psi_forward, double psi_backward);
  void zeta_timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward);
  void jacobis_method_2d_bounded(int n, arma::mat zeta);
  void jacobis_method_2d_periodic(int n, arma::mat zeta);
  void evolve_bounded(bool forwardStep);
  void evolve_periodic(bool forwardStep);
};

#endif
