#include "../include/planet.hpp"

using namespace arma;
using namespace std;

planet::planet(string name_in, double M, vec r, vec v)
{
    name = name_in;
    mass = M;
    c2 = 3992754938.0625; // Squared speed of light in vaccuum [AU/year]
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
    c2 = 3992754938.0625; // Squared speed of light in vaccuum [AU/year]
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

// vec planet::gravitationalForce(const planet& otherPlanet, double Gconst)
// {
//   vec r = distance(otherPlanet);
//   return -Gconst * mass * otherPlanet.mass / pow(norm(r), exp_plus) * r;
// }

vec planet::gravitationalForce(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  double l = specificAngularMomentum(otherPlanet);
  double relcorr = 1 + 3 * pow(norm(l),2) / (pow(norm(r),2) * c2);
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

// double planet::specificAngularMomentum()
// {
//   vec L;
//   if (position.n_elem <= 2){
//     vec r = {position(0), position(1), 0};
//     vec v = {velocity(0), velocity(1), 0};
//     L = cross(r, v);
//   }
//   else{
//     L = cross(position, velocity);
//   }
//   return norm(L);
// }

double planet::kineticEnergy()
{
  return 0.5 * mass * norm(velocity) * norm(velocity);
}

double planet::potentialEnergy(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  return -Gconst * mass * otherPlanet.mass / (exp_minus*pow(norm(r), exp_minus));

}
