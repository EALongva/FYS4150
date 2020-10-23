// testing the functionality of the solar system (ss) class

#include "include/ss.hpp"
#include "include/planet.hpp"
#include "include/utils.hpp"

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
  double C = 365.0;

  // SUN
  std::string name_sun;
  name_sun = "Sun";
  double M_sun;
  double M_sun0;
  M_sun = 2*std::pow(10,30);
  M_sun0 = 1.0;

  arma::vec sun_pos;
  arma::vec sun_vel;
  sun_pos = {0.0,0.0,0.0};
  sun_vel = {0.0,0.0,0.0};


  // EARTH
  std::string name_earth;
  name_earth = "Earth";
  double M_earth;
  M_earth = 6*std::pow(10,24)/M_sun;

  arma::vec earth_pos;
  arma::vec earth_vel;
  earth_pos = {-1.663457613546699*std::pow(10,-1),  9.691203921746886*std::pow(10,-1), -4.125583172010008*std::pow(10,-5)};
  earth_vel = {-1.723919408870981*std::pow(10,-2)*C, -2.981520896064708*std::pow(10,-3)*C,  4.254600200473125*std::pow(10,-7)*C};

  //planet earth(name_earth, M_earth, earth_pos, earth_vel);
  //planet sun(name_sun, M_sun0, sun_pos, sun_vel);

  //earth.print_name();

  //ss solsys;
  //solsys.add(earth);
  //solsys.add(sun);

  //solsys.planet_names();
  //solsys.test(0);


  double T = 1.0;
  int N = std::pow(10.0,N_in);

  ss solsys = import_solarsys();
  solsys.VVerlet(T,N);

  //solsys.VVerlet(T, N);
  //solsys.Xevo[0].print();

  std::cout << "initial pos earth" << std::endl;
  solsys.Xevo[1].col(0).print();
  std::cout << "norm  " << arma::norm(solsys.Xevo[1].col(0)) << std::endl;
  std::cout << std::endl;
  std::cout << "final pos earth" << std::endl;
  solsys.Xevo[1].col(N-1).print(); // this should be implemented as some sort of test
  std::cout << "norm  " << arma::norm(solsys.Xevo[1].col(N-1)) << std::endl;
  std::cout << std::endl;
  std::cout << "initial pos sun" << std::endl;
  solsys.Xevo[0].col(N-1).print();
  std::cout << "norm  " << arma::norm(solsys.Xevo[0].col(N-1)) << std::endl;
  std::cout << std::endl;

  // attempt to save data using arma .save(filename,file_type), csv_ascii for csv format

  for (int i; i < 10; i++){

    number = std::to_string()

    solsys.Xevo[i].save("data/SolarSystem/" + solsys.all_planets[i].name + ".csv")

  }

  solsys.Xevo[0].save("data/EarthSun/earth.csv", arma::csv_ascii);
  //solsys.Xevo[1].save("data/EarthSun/sun.csv", arma::csv_ascii);


  return 0;
}
