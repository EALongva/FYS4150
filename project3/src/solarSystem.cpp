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
  totalMass = 0.0;
  totalKinetic = 0.0;
  totalPotential = 0.0;
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

solarSystem::velocityVerlet(double finalTime, double dt)
{
  // Finding the spatial dimenions from the properties of the planets
  planet& firstplanet = allPlanets[0];
  int n = (int) sizeof(firstplanet.position); // dimension n

  // Finding position of the centre of mass
  vec R (n, fill::zeros);
  for planet in allPlanets { R += planet.mass * planet.position; }
  R = R/totalMass;

  // Using the centre of mass as origin and updating initial positions
  for planet in allPlanets{ planet.position -= R; }

  // Updating the initial velocity of the sun to keep the centre of mass fixed
  vec v0Sun (n, fill::zeros);
  for planet in allPlanets { v0Sun -= planet.mass * planet.velocity; }
  firstplanet.velocity = v0Sun;

  
}
