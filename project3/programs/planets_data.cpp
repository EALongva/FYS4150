// testing the functionality of the solar system (ss) class

#include "include/ss.hpp"
#include "include/planet.hpp"

int main(){


  double AU = 149597870.7;
  double year = 365*24*60*60;
  double MassEarth = 5.97237*std::pow(10.,24);
  double conv = year/AU;
  double C = 365.0;

  // All data credit to NASA, all values with sun at the center at time 00:00 01-01-20207
  // The distances are in units of AU, and velocities in units of AU/day, multiply with C to convert to AU/yr
  // All masses originally in KG (obtained from project description),
  // but dividing by the sun mass so that all masses are in units of the sun mass

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
  earth_vel = {-1.723919408870981*std::pow(10,-2)*C, -2.981520896064708*std::pow(10,-3)*C,  4.254600200473125*std::pow(10,-7)*C}; //<- multiply by C=365


  // MERCURY
  std::string name_mercury;
  name_mercury = "Mercury";
  double M_mercury;
  M_mercury = 3.3*std::pow(10,23)/M_sun;

  arma::vec mercury_pos;
  arma::vec mercury_vel;
  mercury_pos = {-6.333487572394930*std::pow(10,-2), -4.608453269808703*std::pow(10,-1), -3.184761165078634*std::pow(10,-2)};
  mercury_vel = {2.222816779156590*std::pow(10,-2)*C, -2.399853089908365*std::pow(10,-3)*C, -2.235205883702246*std::pow(10,-3)*C};


  // VENUS
  std::string name_venus;
  name_venus = "Venus";
  double M_venus;
  M_venus = 4.9*std::pow(10,24)/M_sun

  arma::vec venus_pos;
  arma::vec venus_vel;
  venus_pos = {7.232002999670082*std::pow(10,-1),  5.254837930328842*std::pow(10,-2), -4.101276234678778*std::pow(10,-2)};
  venus_vel = {-1.547265569012768*std::pow(10,-3)*C,  2.008132769363285*std::pow(10,-2)*C,  3.648328693612742*std::pow(10,-4)*C};


  // MARS
  std::string name_mars;
  name_mars = "Mars";
  double M_mars;
  M_mars = 6.6*std::pow(10,23)/M_sun;

  arma::vec mars_pos;
  arma::vec mars_vel;
  mars_pos = {-1.320107604952232, -8.857574644771996*std::pow(10,-1),  1.382907562643052*std::pow(10,-2)};
  mars_vel = {8.320854741090488*std::pow(10,-3)*C, -1.042277973806052*std::pow(10,-2)*C, -4.225617660650718*std::pow(10,-4)*C};


  // JUPITER
  std::string name_jupiter;
  name_jupiter = "Jupiter";
  double M_jupiter;
  M_jupiter = 1.9*std::pow(10,27)/M_sun;

  arma::vec jupiter_pos;
  arma::vec jupiter_vel;
  jupiter_pos = {5.261470562232079*std::pow(10,-1), -5.201022508399864,  9.830503793253315*std::pow(10,-3)};
  jupiter_vel = {7.423674451944973*std::pow(10,-3)*C,  1.116865602956755*std::pow(10,-3)*C, -1.707572410053752*std::pow(10,-4)*C};


  // SATURN
  std::string name_saturn;
  name_saturn = "Saturn";
  double M_saturn;
  M_saturn = 5.5*std::pow(10,26)/M_sun;

  arma::vec saturn_pos;
  arma::vec saturn_vel;
  saturn_pos = {3.797244866040094, -9.288096505203059,  1.033317176571431*std::pow(10,-2)};
  saturn_vel = {4.863508513301117*std::pow(10,-3)*C,  2.097206392769673*std::pow(10,-3)*C, -2.303054782838117*std::pow(10,-4)*C};


  // URANUS
  std::string name_uranus;
  name_uranus = "Uranus";
  double M_uranus;
  M_uranus = 8.8*std::pow(10,25)/M_sun;

  arma::vec uranus_pos;
  arma::vec uranus_vel;
  uranus_pos = {1.622549696932565*std::pow(10,1),  1.137903567010707*std::pow(10,1), -1.678873902475392*std::pow(10,-1)};
  uranus_vel = {-2.280097293712969*std::pow(10,-3)*C,  3.037887411682148*std::pow(10,-3)*C,  4.074765947142335*std::pow(10,-5)*C};


  // NEPTUNE
  std::string name_neptune;
  name_neptune = "Neptune";
  double M_neptune;
  M_neptune = 1.03*std::pow(10,26)/M_sun;

  arma::vec neptune_pos;
  arma::vec neptune_vel;
  neptune_pos = {2.924281358757302*std::pow(10,1), -6.367188613676754, -5.428996901354922*std::pow(10,-1)};
  neptune_vel = {6.545263007224224*std::pow(10,-4)*C,  3.087626414835625*std::pow(10,-3)*C, -7.892744895737219*std::pow(10,-5)*C};

  // PLUTO
  std::string name_pluto;
  name_pluto = "Pluto";
  double M_pluto;
  M_pluto = 1.31*std::pow(10,22)/M_sun;

  arma::vec pluto_pos;
  arma::vec pluto_vel;
  pluto_pos = {1.297671019650969*std::pow(10,1), -3.137017180006011*std::pow(10,1), -3.965604661998776*std::pow(10,-1)};
  pluto_vel = {2.976602060639467*std::pow(10,-3)*C,  5.376699780750989*std::pow(10,-4)*C, -9.039013018324525*std::pow(10,-4)*C};


  //planet earth(name_earth, M_earth, earth_pos, earth_vel);
  //planet sun(name_sun, M_sun0, sun_pos, sun_vel);


  return 0;
}
