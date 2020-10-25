#ifndef PLANET_H
#define PLANET_H

#include <iostream>
#include <cmath>
#include <armadillo>

class planet
{
public:
  // Planet Properties
  std::string name;
  double mass;
  double c2;
  arma::vec position;
  arma::vec velocity;
  double potential;
  double kinetic;
  double force;
  double exponent;
  double exp_plus;
  double exp_minus;

  // Planet Initializer
  planet(std::string name_in, double M, arma::vec r, arma::vec v);
  planet(std::string name_in, double M, arma::vec r, arma::vec v, double beta);

  // Functions
  void print_name();
  arma::vec distance(const planet& otherPlanet);
  arma::vec gravitationalForce(const planet& otherPlanet, double Gconst);
  arma::vec gravitationalForce_relcorr(const planet& otherPlanet, double Gconst);
  arma::vec acceleration(const planet& otherPlanet, double Gconst);
  double specificAngularMomentum(const planet& otherPlanet);
  double specificAngularMomentum2();
  double kineticEnergy();
  double potentialEnergy(const planet& otherPlanet, double Gconst);
};

#endif
