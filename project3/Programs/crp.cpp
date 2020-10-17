// testing the functionality of the solar system (ss) class

#include "include/ss.hpp"
#include "include/planet.hpp"

int main(int argc, char* argv[]){

  int N_in;
  N_in = std::atoi(argv[1]);

  /*
  double fx;
  double dfx;
  fx = 5;
  dfx = 2;

  ss test(fx);
  test.print();
  test.change_var(dfx);
  test.print();
  */


  // testing planet class (earth)

  double AU = 149597870.7;
  double year = 365*24*60*60;
  double MassEarth = 5.97237*std::pow(10.,24);
  double conv = year/AU;

  arma::vec r_0;
  arma::vec v_0;
  double vy_0;
  double M_earth;
  double M_sun;
  std::string name_earth;
  std::string name_mars;
  std::string name_sun;

  name_mars = "Mars";
  name_earth = "Earth";
  name_sun = "Sun";
  M_earth = 1.0/332933.;
  M_sun = 1.0;
  vy_0 = -30.29*conv;
  r_0 = {-0.98329, 0.0};
  v_0 = {0.0, vy_0};
  arma::vec sunpos = {0.0, 0.0};
  arma::vec sunvel = {0.0, 0.0};

  planet earth(name_earth, M_earth, r_0, v_0);
  planet mars(name_mars, M_earth, r_0, v_0);
  planet sun(name_sun, M_sun, sunpos, sunvel);

  earth.print_name();

  ss solsys;

  solsys.add(earth);
  solsys.add(sun);

  //solsys.planet_names();
  //solsys.test(0);


  double T = 1.0;
  int N = std::pow(10.0,N_in);
  solsys.VVerlet(T, N);
  //solsys.Xevo[0].print();
  solsys.Xevo[0].col(0).print();
  solsys.Xevo[0].col(N-1).print(); // this should be implemented as some sort of test
  solsys.Xevo[1].col(N-1).print();

  return 0;
}
