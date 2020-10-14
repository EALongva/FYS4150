// testing the functionality of the solar system (ss) class

#include "include/ss.hpp"
#include "include/planet.hpp"

int main(){

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

  arma::vec r_0;
  arma::vec v_0;
  double vx_0;
  double M_earth;
  std::string name_earth;

  name_earth = "Earth";
  M_earth = 1.0;
  vx_0 = M_PI;
  r_0 = {0.0, 1.0};
  v_0 = {vx_0, 0.0};

  planet earth(name_earth, M_earth, r_0, v_0);

  earth.print_name();




  return 0;
}
