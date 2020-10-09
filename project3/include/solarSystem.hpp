#include "planet.h"
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
  std::map<string,int> planetIndices;
  double radius;
  double totalMass;
  double totalKinetic;
  double totalPotential;
s
  // Solar System Initializer
  solarSystem(double rad);

  // Functions
  void add_planet(const planet& newPlanet);
  planet get_planet(std::string planetName);
}
