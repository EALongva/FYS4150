#ifndef ROSSBY_H
#define ROSSBY_H

#include <iostream>
#include <cmath>
#include <armadillo>

class rossby
{
public:
  // Properties
  double endpos;
  double endtime;
  double deltax;
  double deltat;
  int xdim;
  int tdim;
  arma::mat Psi;
  arma::mat Zeta;

  // Initializer
  rossby(double dx, double dt, double endtime);
  // rossby(double u_in, double v_in);
  // rossby(double u_in, double v_in, arma::vec psi_in);


  // Functions
  //void function();
  void initialize_wave(bool sineWave, double sigma, double x0);
  void zeta_timestep_forward(double &zeta_forward, double psi_forward, double psi_backward);
  void zeta_timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward);
  arma::vec precalculate_offdiag();
  arma::vec forward_substitution(arma::vec zeta, arma::vec c_new);
  void evolve_bounded(bool forwardStep);
};

#endif
