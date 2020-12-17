#include <iostream>
#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "time.h"
#include "include/rossby.hpp"

using namespace arma;
using namespace std;

int main(int argc, char *argv[]){
  string psiname = "../results/data/psi_2d_bounded";

  // Takes size of grid space, size of time step and length of time period as input
  double deltapos = 0.01;
  double deltatime = 0.01;
  double endtime = 2*8*M_PI*M_PI;
  double sigma = 0.1; double x0 = 0.5; double y0 = 0.5;

  // analytical phase speed c
  double L = 1;
  double beta = 1;
  double n = 1;

  double c_p = -(beta*L*L) / (4*n*n*M_PI*M_PI);
  double c_b = -(beta*L*L) / (2*n*n*M_PI*M_PI);

  // analytic stream functions
  int nx;
  nx = (int) 1.0/deltapos;
  int nt;
  nt = (int) endtime/deltatime;
  double A = 1;

  vec psiPer = vec(nx,fill::zeros);
  vec psiPer0 = vec(nx,fill::zeros);
  vec psiPerHalf = vec(nx,fill::zeros);

  double k = 4*M_PI;
  double omega = 1/k;

  for (int i=0;i<nx;i++){
    double x = (i+1)*deltapos;
    psiPer(i) = A*sin(k*x - omega*(endtime));
  }

  for (int i=0;i<nx;i++){
    double x = (i+1)*deltapos;
    psiPer0(i) = A*sin(k*x);
  }

  for (int i=0;i<nx;i++){
    double x = (i+1)*deltapos;
    psiPerHalf(i) = A*sin(k*x - omega*(endtime/2.0));
  }


  //cout << "analytical" << endl;
  //psiPer.print();


  //vec psiBou =


  // calculate phase speed numerically 1D0
  // init sinewave

  bool sineWave = true;
  rossby RB(deltapos, deltatime, endtime);
  RB.initialize_wave(sineWave, sigma, x0);
  //RB.Psi.col(0).print();

  bool forwardStep = false;
  RB.evolve_periodic(forwardStep);

  //cout << "numerical" << endl;
  //RB.Psi.col(nx).print();
  vec num_psiPer = RB.Psi.col(nt-1);
  vec num_psiPer0 = RB.Psi.col(0);
  vec num_psiPerHalf = RB.Psi.col(nt-((int)nt/2));


  vec relError_psiPer = vec(nx,fill::zeros);
  for (int j=0;j<nx;j++){
  relError_psiPer(j) = fabs((num_psiPer(j) - psiPer(j)));
  }

  //cout << "relative error" << endl;

  psiPer.save("../results/data/test_analytical_periodic.csv", csv_ascii);
  num_psiPer.save("../results/data/test_numerical_periodic.csv", csv_ascii);
  relError_psiPer.save("../results/data/test_relerror_periodic.csv", csv_ascii);

  psiPer0.save("../results/data/test_analytical_periodic0.csv", csv_ascii);
  num_psiPer0.save("../results/data/test_numerical_periodic0.csv", csv_ascii);
  psiPerHalf.save("../results/data/test_analytical_periodicHalf.csv", csv_ascii);
  num_psiPerHalf.save("../results/data/test_numerical_periodicHalf.csv", csv_ascii);


  return 0;

}
