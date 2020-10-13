#include "../include/solarSystem.hpp"
#include "../include/planet.hpp"
#include <iostream>
#include <cmath>
#include "time.h"
#include <typeinfo>

using namespace arma;
using namespace std;

solarSystem::solarSystem(double rad, double Gconst)
{
  totalPlanets = 0;
  radius = rad;
  totalMass = 0.0;
  totalKinetic = 0.0;
  totalPotential = 0.0;
}

void solarSystem::add_planet(const planet& newPlanet)
{
  totalPlanets += 1;
  totalMass += newPlanet.mass;
  allPlanets.push_back(newPlanet);
  planetIndices.insert(pair<string,int>(newPlanet.name,(int)totalPlanets-1));
}

planet solarSystem::get_planet(string planetName)
{
  int i = planetIndices.at(planetName);
  return allPlanets[i];
}

void solarSystem::velocityVerlet(double finalTime, double dt)
{
  // Finding the spatial dimenions from the properties of the planets
  planet Sun = get_planet("Sun");
  int n = (int) sizeof(Sun.position); // dimension n

  // Finding position of the centre of mass
  vec R (n, fill::zeros);
  for (planet planet:allPlanets){ R += planet.mass * planet.position;}
  R = R/totalMass;

  // Using the centre of mass as origin and updating initial positions
  for (planet planet:allPlanets){ planet.position -= R; }

  // Updating the initial velocity of the sun to keep the centre of mass fixed
  vec v0Sun (n, fill::zeros);
  for (planet planet:allPlanets){ v0Sun -= planet.mass * planet.velocity; }
  Sun.velocity = v0Sun;

  // Initializing force and acceleration vectors
  double time = 0.0;
  vec F(n,fill::zeros); vec F0(n,fill::zeros);
  vec a(n,fill::zeros);
  vec a_new(n,fill::zeros);
  double dt_sqrd = 0.5*dt*dt;
  double dt_half = 0.5*dt;



  // Looping over all planets in time
  while (time < finalTime) {
    time += dt;
    for (planet planet:allPlanets){
      F = F0; // Setting the gravity forces equal to zero
      for (planet otherPlanet:allPlanets){
        F += planet::gravitationalForce(otherPlanet, Gconst);
      }

      a = F/planet.mass;
      planet.position += dt*planet.velocity + dt_sqrd*a;

      F = F0; // Setting the gravity forces equal to zero
      for (planet otherPlanet:allPlanets){
        F += planet::gravitationalForce(otherPlanet, Gconst);
      }

      a_new = F/planet.mass;
      planet.velocity += dt_half*(a_new + a);
    }
  }
}
