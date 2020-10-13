#include "../include/planet.hpp"
#include <vector>
#include <map>
#include <armadillo>

class solarSystem
{
public:
  friend class planet;

  // Solar System Properties
  int totalPlanets;
  std::vector<planet> allPlanets;
  std::map<std::string,int> planetIndices;
  double radius;
  double totalMass;
  double totalKinetic;
  double totalPotential;

  // Solar System Initializer
  solarSystem(double rad, double Gconst);

  // Functions
  void add_planet(const planet& newPlanet);
  planet get_planet(std::string planetName);
  void velocityVerlet(double finalTime, double dt);
};
