// class for project 3

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "time.h"
#include "../include/planet.hpp"
#include "../include/ss.hpp"

using namespace arma;
using namespace std;

ss::ss()
{
  //G = 4*M_PI*M_PI;
  //G = 0.000011;
  double AU = 149597870.7;
  double year = 365*24*60*60;
  double MassEarth = 5.97237*std::pow(10.,24);
  double conv = year/AU;
  double conv2 = MassEarth*conv*conv/AU;

  double G_kiloSI = 6.67430*std::pow(10.,-20);

  //G = G_kiloSI*conv2;
  G = 4*M_PI*M_PI;

  total_planets = 0;
}

void ss::add(planet newplanet)
{
  total_planets += 1;
  all_planets.push_back(newplanet);
}

void ss::planet_names()
{
  for (int i=0; i < total_planets; i++){

    cout << all_planets[i].name << ", " << i;

  }

  cout << endl;

}


void ss::euler(double T, int N)
{

  std::cout << G << std::endl;

  double dt = T / ((double) N);

  int dim = all_planets[0].position.n_elem;

  arma::mat M( dim, N, arma::fill::zeros);

  for (int nr = 0; nr < total_planets; nr++){

    M.col(0) = all_planets[nr].position;
    Xevo.push_back(M);
    M.col(0) = all_planets[nr].velocity;
    Vevo.push_back(M);

  }

  arma::mat A( dim, total_planets, arma::fill::zeros);

  for (int i = 0; i < N-1; i++){

    for (int nr1 = 0; nr1 < total_planets; nr1++){
      for (int nr2 = 0; nr2 < total_planets; nr2++){

        if ( nr1 != nr2 ){
          A.col(nr1) += all_planets[nr1].acceleration(all_planets[nr2], G);
        }

      }
    }

    for (int j = 0; j < total_planets; j++){

      Xevo[j].col(i+1) = Xevo[j].col(i) + dt*Vevo[j].col(i);
      Vevo[j].col(i+1) = dt*A.col(j);

      all_planets[j].position = Xevo[j].col(i+1);
      // reset acceleration matrix A
      A.col(j) = A.col(j)*0;

    }

  }

}
//*/

void ss::eulerStep(int I, int J, int i, double dt)
{

  Xevo[I].col(i+1) = Xevo[I].col(i) + dt*Vevo[I].col(i);
  Vevo[I].col(i+1) = dt*all_planets[I].acceleration(all_planets[J], G);

}

void ss::VVerlet(double T, int N)
{

  std::cout << G << std::endl;

  double dt = T / ((double) N);

  int dim = all_planets[0].position.n_elem;

  arma::mat M( dim, N, arma::fill::zeros);

  for (int nr = 0; nr < total_planets; nr++){

    M.col(0) = all_planets[nr].position;
    Xevo.push_back(M);
    M.col(0) = all_planets[nr].velocity;
    Vevo.push_back(M);

  }

  arma::mat A1( dim, total_planets, arma::fill::zeros);
  arma::mat A2( dim, total_planets, arma::fill::zeros);

  for (int i = 0; i < N-1; i++){

    for (int nr1 = 0; nr1 < total_planets; nr1++){
      for (int nr2 = 0; nr2 < total_planets; nr2++){

        if ( nr1 != nr2 ){
          A1.col(nr1) += all_planets[nr1].acceleration(all_planets[nr2], G);
        }

      }
    }

    for (int j = 0; j < total_planets; j++){

      Xevo[j].col(i+1) = Xevo[j].col(i) + dt*Vevo[j].col(i) + dt*dt*0.5*A1.col(j);
      //Vevo[j].col(i+1) = dt*A1.col(j);

      all_planets[j].position = Xevo[j].col(i+1);

    }

    for (int nr1 = 0; nr1 < total_planets; nr1++){
      for (int nr2 = 0; nr2 < total_planets; nr2++){

        if ( nr1 != nr2 ){
          A2.col(nr1) += all_planets[nr1].acceleration(all_planets[nr2], G);
        }

      }
    }

    for (int j = 0; j < total_planets; j++){

      Vevo[j].col(i+1) = Vevo[j].col(i) + dt*0.5*(A2.col(j) + A1.col(j));
      // reset acceleration matrix A
      A1.col(j) = A1.col(j)*0;
      A2.col(j) = A2.col(j)*0;


    }


  }

}


void ss::test(int nr)
{

  //std::cout << all_planets[nr].position.n_elem << std::endl;

  int dim = all_planets[0].position.n_elem;

  arma::mat M( dim, 10, arma::fill::zeros);

  M.col(0) = all_planets[nr].position;
  M.print();

  std::cout << all_planets[nr].position << std::endl;
  all_planets[nr].position = {3.0, 4.0};
  std::cout << all_planets[nr].position << std::endl;

  double hello = (2*3)/2;
  std::cout << hello << std::endl;


}



//
