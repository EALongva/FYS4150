#include "../include/solarSystem.hpp"
#include "../include/planet.hpp"
#include <iostream>
#include <cmath>
#include "time.h"
#include <typeinfo>

using namespace arma;
using namespace std;

solarSystem::solarSystem(int dim, double G, double rad)
{
  dimension = dim;
  totalPlanets = 0;
  radius = rad;
  Gconst = G;
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

void solarSystem::FixOriginCentreOfMass()
{
  planet Sun = get_planet("Sun");

  // Updating the initial velocity of the Sun to make the total momentum zero
  // (to keep the centre of mass fixed)
  vec v0Sun (dimension, fill::zeros);
  for (planet planet:allPlanets){ v0Sun -= planet.mass * planet.velocity; }
  Sun.velocity = v0Sun/Sun.mass;

  // Finding the position of the centre of mass
  vec R (dimension, fill::zeros);
  for (planet planet:allPlanets){ R += planet.mass * planet.position;}
  R = R/totalMass;

  // Using the centre of mass as origin and updating initial positions
  for (planet planet:allPlanets){ planet.position -= R; }
}

vec solarSystem::acceleration(int i)
{
  planet planet = allPlanets[i];
  vec a(dimension,fill::zeros);
  vec F(dimension,fill::zeros);
  vec zerovec(dimension, fill::zeros);
  for (int j = 0; j < totalPlanets; j++){
    if (i == j){
      F += zerovec;
    }
    else{
      F += planet.gravitationalForce(allPlanets[j], Gconst);
    }
  }
  a = F/planet.mass;
  return a;
}

void solarSystem::velocityVerlet(double finalTime, int integrationPoints)
{
  double dt = finalTime/integrationPoints; // size of time step

  mat a(totalPlanets, dimension, fill::zeros);
  mat a_new(dimension, totalPlanets, fill::zeros);
  mat r(dimension, totalPlanets, fill::zeros);
  mat v(dimension, totalPlanets, fill::zeros);
  double dt_pos = 0.5*dt*dt;
  double dt_vel = 0.5*dt;

  for (int i = 0; i < totalPlanets; i++){
    planet planet = allPlanets[i];
    a.col(i) = acceleration(i);
    planet.position += dt*planet.velocity + dt_sqrd*a.col(i);
  }
  for (int i = 0; i < totalPlanets, i++){
    
  }
}


/*
void solarSystem::velocityVerlet(double finalTime, int integrationPoints)
{
  double dt = finalTime/integrationPoints; // size of time step

  // Initializing force and acceleration vectors
  double time = 0.0;
  mat F(N,n,fill::zeros);
  //mat a(n,fill::zeros);
  //mat a_new(n,fill::zeros);
  double dt_sqrd = 0.5*dt*dt;
  double dt_half = 0.5*dt;

  for (int i = 0; i < N; i++){
    planet planet = allPlanets[i];
    F.row(i).fill(0.0); // Setting the gravity forces equal to zero
    for (int j = 0; j < N; j++){
      planet otherPlanet = allPlanets[j];
      F.row(i) += planet::gravitationalForce(otherPlanet, Gconst);
    }








  // Looping over all planets in time
  while (time < finalTime) {
    time += dt;
    for (planet planet:allPlanets){
      F.fill(0.0); // Setting the gravity forces equal to zero
      for (planet otherPlanet:allPlanets){
        F += (vec) planet::gravitationalForce(otherPlanet, Gconst);
      }
      a = (vec) F/planet.mass;
      planet.position += dt*planet.velocity + dt_sqrd*a;
      F.fill(0.0); // Setting the gravity forces equal to zero
      for (planet otherPlanet:allPlanets){
        F += planet::gravitationalForce(otherPlanet, Gconst);
      }

      a_new = F/planet.mass;
      planet.velocity += dt_half*(a_new + a);
    }
  }
}
*/
