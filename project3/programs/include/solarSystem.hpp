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
  arma::vec acceleration_relcorr(int i);
  double potentialEnergy(int i);
  void fixOriginCentreOfMass();
  void velocityVerlet(double finalTime, int integrationPoints, std::string outputfilename);
  void perihelionAngle(double finalTime, int integrationPoints, std::string outputfilename);
  void perihelionAngle2(double finalTime, int integrationPoints, std::string outputfilename);
  void perihelionAngle_relcorr(double finalTime, int integrationPoints, std::string outputfilename);
  void perihelionAngle_relcorr2(double finalTime, int integrationPoints, std::string outputfilename);
};
