#include "../include/solarSystem.hpp"
#include "../include/planet.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "time.h"

using namespace arma;
using namespace std;

ofstream ofile;

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

void solarSystem::fixOriginCentreOfMass()
{
  // Finding the position of the centre of mass
  vec R (dimension, fill::zeros);
  for (planet& planet:allPlanets){ R += planet.mass * planet.position;}
  R = R/totalMass;

  // Using the centre of mass as origin and updating initial positions
  for (planet& planet:allPlanets){ planet.position -= R; }

  planet Sun = get_planet("Sun");
  // Updating the initial velocity of the Sun to make the total momentum zero
  // (to keep the centre of mass fixed)
  vec v0Sun (dimension, fill::zeros);
  for (planet& planet:allPlanets){ v0Sun -= planet.mass * planet.velocity; }
  Sun.velocity = v0Sun/Sun.mass;

}

vec solarSystem::acceleration(int i)
{
  planet& planet = allPlanets[i];
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

double solarSystem::potentialEnergy(int i)
{
  planet& planet = allPlanets[i];
  double PE = 0.0;
  for (int j = 0; j < totalPlanets; j++){
    if (i != j){
      PE += planet.potentialEnergy(allPlanets[j], Gconst);
    }
  }
  return PE;
}

void solarSystem::velocityVerlet(double finalTime, int integrationPoints, string outputfilename)
{
  double dt = finalTime/integrationPoints; // Size of time step

  mat A(totalPlanets, dimension, fill::zeros); // Store acceleration for every planet each time step
  mat A_new(totalPlanets, dimension, fill::zeros); // Store new acceleration for every planet each time step
  double PE; // Store potential energy for each step
  vec a(dimension, fill::zeros); // Store output of acceleration function
  double dt_pos = 0.5*dt*dt;
  double dt_vel = 0.5*dt;
  double time = 0.0;

  // Open file for output values
  ofile.open(outputfilename);
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << time << endl;

  // Find the initial acceleration
  for (int i = 0; i < totalPlanets; i++){
    planet& planet = allPlanets[i];
    a = acceleration(i);
    for (int k = 0; k < dimension; k++){
      A_new(i,k) = a(k); // Store as new acceleration
      ofile << planet.position[k]<< endl; // Write out initial planet position
    }
  }

  for (int i = 0; i < totalPlanets; i++){
    planet& planet = allPlanets[i];
    ofile << planet.kineticEnergy() << endl; // Write out initial kinetic energy
    PE = potentialEnergy(i);
    ofile << PE << endl; // Write out initial potential energy
    //ofile << planet.specificAngularMomentum(allPlanets[0]) << endl; // Write out angular momentum
    ofile << planet.specificAngularMomentum() << endl;
  }

  while (time < finalTime){ // Loop over each time step
    time += dt;
    ofile << time << endl;

    for (int i = 0; i < totalPlanets; i++){
      planet& planet = allPlanets[i];
      for (int k = 0; k < dimension; k++){
        A(i,k) = A_new(i,k); // Set current acceleration to new acceleration from previous time step
        planet.position(k) += dt*planet.velocity(k) + dt_pos*A(i,k); // Calculate current position
        ofile << planet.position[k]<< endl; // Write out planet position
      }
    }

    for (int i = 0; i < totalPlanets; i++){
      planet& planet = allPlanets[i];
      PE = potentialEnergy(i);
      a = acceleration(i); // Calculate new acceleration based on current position
      for (int k = 0; k < dimension; k++){
        A_new(i,k) = a(k);
        planet.velocity(k) += dt_vel*(A_new(i,k) + A(i,k)); // Calculate current velocity
      }
      ofile << planet.kineticEnergy() << endl;
      ofile << PE << endl;
      //ofile << planet.specificAngularMomentum(allPlanets[0]) << endl;
      ofile << planet.specificAngularMomentum() << endl;
    }
  }
  ofile.close();
}
