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
}

void planet::print_name()
{
  cout << name << endl;
}

vec planet::distance(const planet& otherPlanet)
{
  return otherPlanet.position - position;
}

vec planet::gravitationalForce(const planet& otherPlanet, double Gconst)
{
  vec r = distance(otherPlanet);
  if (norm(r)!=0){
    return -Gconst * mass * otherPlanet.mass / pow(norm(r),3) * r;
  }
  else {return 0;}
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
  if(norm(r)!=0){
    return -Gconst * mass * otherPlanet.mass / norm(r);
  }
  else {return 0;}
}
