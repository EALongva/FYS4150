#ifndef ROSSBY_2D_H
#define ROSSBY_2D_H

#include <iostream>
#include <cmath>
#include <armadillo>

class rossby
{
public:
  // Properties
  double endposx;
  double endposy;
  double endtime;
  double deltax;
  double deltay;
  double deltat;
  int xdim;
  int ydim;
  int tdim;
  arma::cube Psi;
  arma::cube Zeta;

  // Initializer
  rossby(double dpos, double dt, double endtime);
  // rossby(double u_in, double v_in);
  // rossby(double u_in, double v_in, arma::vec psi_in);


  // Functions
  //void function();
  void initialize_wave(bool sineWave, double sigma, double x0);
  void zeta_timestep_forward(double &zeta_forward, double zeta, double psi_forward, double psi_backward);
  void zeta_timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward);
  void jacobis_method_2d(int n, arma::vec zeta);
  void evolve_bounded(bool forwardStep);
  void evolve_periodic(bool forwardStep);
};

#endif
