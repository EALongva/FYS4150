#include "solarSystem.hpp"
#include "planet.hpp"
#include <iostream>
#include <cmath>
#include "time.h"

using namespace arma;
using namespace std;

solarSystem::solarSystem(double rad)
{
  totalPlanets = 0;
  radius = rad;
  totalMass = 0;
  totalKinetic = 0;
  totalPotential = 0;
}

solarSystem::add_planet(newPlanet)
{
  totalPlanets += 1;
  totalMass += newPlanet.mass;
  allPlanets.push_back(newPlanet);

}

solarSystem::get_planet(planetName)
{
  int i = planetIndices(planetName);
  return allPlanets(i);
}

solarSystem::velocityVerlet(int dimension, double finalTime, double dt)
{
  
}
