#include "../include/planet.hpp"

using namespace arma;
using namespace std;

planet::planet(string name_in, double M, vec r, vec v)
{
    name = name_in;
    mass = M;
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

vec planet::acceleration(const planet& otherPlanet, double Gconst)
{
  vec gForce = gravitationalForce(otherPlanet, Gconst);
  return gForce / mass;
}

double planet::kineticEnergy()
{
  return 0.5 * mass * norm(velocity) * norm(velocity);
}

double planet::potentialEnergy(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  return -Gconst * mass * otherPlanet.mass / exp_minus*pow(norm(r), exp_minus);

}
