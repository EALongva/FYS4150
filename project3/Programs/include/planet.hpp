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
  arma::vec position;
  arma::vec velocity;
  double potential;
  double kinetic;
  double force;

  // Planet Initializer
  planet(std::string name_in, double M, arma::vec r, arma::vec v);

  // Functions
  void print_name();
  arma::vec distance(const planet& otherPlanet);
  arma::vec gravitationalForce(const planet& otherPlanet, double Gconst);
  arma::vec acceleration(const planet& otherPlanet, double Gconst);
  double kineticEnergy();
  double potentialEnergy(const planet& otherPlanet, double Gconst);
};

#endif
