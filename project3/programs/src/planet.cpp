#include "../include/planet.hpp"

using namespace arma;
using namespace std;

planet::planet(string name_in, double M, vec r, vec v)
{
    name = name_in;
    mass = M;
    c2 = 3987291025; // Squared speed of light in vaccuum [AU²/year²]
    position = r;
    velocity = v;
    kinetic = 0.0;
    potential = 0.0;
    force = 0.0;
    exponent = 2;
    exp_plus = exponent+1;
    exp_minus = exponent-1;
}

planet::planet(string name_in, double M, vec r, vec v, double beta)
{
    name = name_in;
    c2 = 3987291025; // Squared speed of light in vaccuum [AU²/year²]
    mass = M;
    position = r;
    velocity = v;
    kinetic = 0.0;
    potential = 0.0;
    force = 0.0;
    exponent = beta;
    exp_plus = exponent+1;
    exp_minus = exponent-1;
}

void planet::print_name()
{
  cout << name << endl;
}

vec planet::distance(const planet& otherPlanet)
{
  return position - otherPlanet.position;
}

vec planet::gravitationalForce(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  return -Gconst * mass * otherPlanet.mass / pow(norm(r), exp_plus) * r;
}

vec planet::gravitationalForce_relcorr(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  //double l = specificAngularMomentum(otherPlanet);
  double l = specificAngularMomentum(otherPlanet);
  double relcorr = 1 + ((3*l*l) / (pow(norm(r),2) * c2));
  return -Gconst * mass * otherPlanet.mass / pow(norm(r), exp_plus) * relcorr * r;
}

vec planet::acceleration(const planet& otherPlanet, double Gconst)
{
  vec gForce = gravitationalForce(otherPlanet, Gconst);
  return gForce / mass;
}

double planet::specificAngularMomentum(const planet& otherPlanet)
{
  vec L;
  vec r = distance(otherPlanet);
  if (r.n_elem <= 2){
    vec rnew = {r(0), r(1), 0};
    vec vnew = {velocity(0), velocity(1), 0};
    L = cross(rnew, vnew);
  }
  else{
    L = cross(r, velocity);
  }
  return norm(L);
}

double planet::kineticEnergy()
{
  return 0.5 * mass * norm(velocity) * norm(velocity);
}

double planet::potentialEnergy(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  return -Gconst * mass * otherPlanet.mass / (exp_minus*pow(norm(r), exp_minus));

}
