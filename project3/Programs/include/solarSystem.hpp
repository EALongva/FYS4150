#include "../include/planet.hpp"
#include <vector>
#include <map>
#include <armadillo>

class solarSystem
{
public:
  friend class planet;

  // Solar System Properties
  int dimension;
  int totalPlanets;
  std::vector<planet> allPlanets;
  std::map<std::string,int> planetIndices;
  double radius;
  double totalMass;
  double totalKinetic;
  double totalPotential;
  double Gconst;

  // Solar System Initializer
  solarSystem(int dim, double G, double rad);

  // Functions
  void add_planet(const planet& newPlanet);
  planet get_planet(std::string planetName);
  arma::vec acceleration(int i);
  double potentialEnergy(int i);
  void fixOriginCentreOfMass();
  void velocityVerlet(double finalTime, int integrationPoints, std::string outputfilename);
};
