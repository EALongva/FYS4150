
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include "time.h"
using namespace  std;
using namespace arma;
// output file
ofstream ofile;

// inline function for PeriodicBoundary boundary conditions
inline int PeriodicBoundary(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void InitializeLattice(int, mat &, double&, double&);
// The metropolis algorithm including the loop over Monte Carlo cycles
void MetropolisSampling(int, int, double, vec &);
// prints to file the results of the calculations
void WriteResultstoFile(int, int, double, vec);

// Main program begins here

int main(int argc, char* argv[])
{
  string filename, Orientation;
  int NSpins, MonteCarloCycles;
  double InitialTemp, FinalTemp, TempStep;
  if (argc <= 5) {
    cout << "Bad Usage: " << argv[0] <<
      " read output file, Number of spins, MC cycles, temperature and orientation ('random' or 'ordered')" << endl;
    exit(1);
  }

  string outpath = "results/data/";
  filename = outpath+argv[1];
  NSpins = atoi(argv[2]);
  MonteCarloCycles = atoi(argv[3]);
  Temperature = atof(argv[4]);
  Orientation = argv[5]


  // Declare new file name and add temperature and orientation to file name, only master node opens file
  string fileout = filename;
  string argument = to_string(Temperature);
  fileout.append(argument+Orientation);
  ofile.open(fileout);

  // Start Monte Carlo sampling by looping over the selected Temperatures
  double  TimeStart, TimeEnd, TotalTime;
  clock_t TimeStart = clock();

  vec AllEnergies = zeros<mat>(MonteCarloCycles);
  vec AllMagnetisations = zeros<mat>(MonteCarloCycles);

  // Start Monte Carlo computation and get expectation values
  MetropolisSampling(NSpins, MonteCarloCycles, Temperature, AllEnergies, AllMagnetisations);

  // Find total average
  WriteResultstoFile(NSpins, AllEnergies, AllMagnetisations);
  ofile.close();  // close output file
  clock_t TimeEnd = clock();
  TotalTime = (TimeStart-TimeEnd)/CLOCKS_PER_SEC);
  cout << "Time = " <<  TotalTime  << endl;
  // End MPI
  return 0;
}


// The Monte Carlo part with the Metropolis algo with sweeps over the lattice
void MetropolisSampling(int NSpins, int MonteCarloCycles, double Temperature, vec &AllEnergies, vec &AllMagnetisations)
{
  // Initialize the seed and call the Mersienne algo
  std::random_device rd;
  std::mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in [[0, 1]
  std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  // Initialize the lattice spin values
  mat SpinMatrix = zeros<mat>(NSpins,NSpins);
  //    initialize energy and magnetization
  double Energy = 0.;     double MagneticMoment = 0.;
  // initialize array for expectation values
  InitializeLattice(NSpins, SpinMatrix, Energy, MagneticMoment, Orientation);
  // setup array for possible energy changes
  vec EnergyDifference = zeros<mat>(17);
  for( int de =-8; de <= 8; de+=4) EnergyDifference(de+8) = exp(-de/Temperature);
  // Start Monte Carlo experiments
  int AllSpins = NSpins*NSpins;
  for (int cycles = 1; cycles <= MonteCarloCycles; cycles++){
    // The sweep over the lattice, looping over all spin sites
    for(int Spins =0; Spins < AllSpins; Spins++) {
      int ix = (int) (RandomNumberGenerator(gen)*NSpins);
      int iy = (int) (RandomNumberGenerator(gen)*NSpins);
      int deltaE =  2*SpinMatrix(ix,iy)*
	(SpinMatrix(ix,PeriodicBoundary(iy,NSpins,-1))+
	 SpinMatrix(PeriodicBoundary(ix,NSpins,-1),iy) +
	 SpinMatrix(ix,PeriodicBoundary(iy,NSpins,1)) +
	 SpinMatrix(PeriodicBoundary(ix,NSpins,1),iy));
      if ( RandomNumberGenerator(gen) <= EnergyDifference(deltaE+8) ) {
	SpinMatrix(ix,iy) *= -1.0;  // flip one spin and accept new spin config
	MagneticMoment += 2.0*SpinMatrix(ix,iy);
	Energy += (double) deltaE;
      }
  }
  AllEnergies(cycles) = Energy;
  AllMagnetisations(cycles) = fabs(MagneticMoment);
} // end of Metropolis sampling over spins

// function to initialise energy, spin matrix and magnetization
void InitializeLattice(int NSpins, mat &SpinMatrix,  double& Energy, double& MagneticMoment, string Orientiation)
{
  if (Orientation == "random"){
    // Initialize the seed and call the Mersienne algo
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Set up the uniform distribution for x \in [[0, 1]
    std::uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    // setup spin matrix and initial magnetization using random orientation
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
        if ( RandomNumberGenerator(gen) <= 0.5 ) {
          SpinMatrix(x,y) = -1.0; // spin orientation for the ground state
          MagneticMoment +=  (double) SpinMatrix(x,y);
        }
        else {
          SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
          MagneticMoment +=  (double) SpinMatrix(x,y);
        }
      }
    }
  }
  else{
    // setup spin matrix and initial magnetization using cold start, all spins pointing up or down
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
        SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
        MagneticMoment +=  (double) SpinMatrix(x,y);
      }
    }
  }
    // setup spin matrix and initial magnetization using cold start, all spins pointing up or down
    for(int x =0; x < NSpins; x++) {
      for (int y= 0; y < NSpins; y++){
        SpinMatrix(x,y) = 1.0; // spin orientation for the ground state
        MagneticMoment +=  (double) SpinMatrix(x,y);
      }
    }
  }
  // setup initial energy
  for(int x =0; x < NSpins; x++) {
    for (int y= 0; y < NSpins; y++){
      Energy -=  (double) SpinMatrix(x,y)*
	(SpinMatrix(PeriodicBoundary(x,NSpins,-1),y) +
	 SpinMatrix(x,PeriodicBoundary(y,NSpins,-1)));
    }
  }
}// end function initialize



void WriteResultstoFile(int NSpins, vec AllEnergies, vec AllMagnetisations)
{
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << AllEnergies
  ofile << setw(15) << setprecision(8) << AllMagnetisations << endl;
} // end output function