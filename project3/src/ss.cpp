// class for project 3

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "time.h"
#include "../include/planet.hpp"
#include "../include/ss.hpp"

ss::ss()
{
  G = 4*M_PI*M_PI;
  total_planets = 0;
}

void ss::add(planet newplanet)
{
  total_planets ++;
  all_planets.push_back(newplanet);
}

void ss::planet_names()
{
  for (int i=0; i < )
}




//
