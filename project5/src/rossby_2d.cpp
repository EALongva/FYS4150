#include "../include/rossby_2d.hpp"

using namespace arma;
using namespace std;

inline int periodic(int i, int limit, int add){
  return (i+limit+add) % (limit);}


//  Initializers

rossby::rossby(double dpos, double dt, double tfinal)
{
  endposx = 1.0;
  endposy = 1.0;
  endtime = tfinal;
  xdim = (int) endposx/dpos;
  ydim = (int) endposy/dpos;
  tdim = (int) endtime/dt;
  deltax = dpos;
  deltay = dpos;
  deltat = dt;
  Psi = cube(xdim, ydim, tdim, fill::zeros);
  Zeta = cube(xdim, ydim, tdim, fill::zeros);
  cout << "rossby class object initialized successfully" << endl;
}
//
// rossby::rossby(double u_in, double v_in)
// {
//   u = u_in;
//   v = v_in;
//   std::cout << "rossby class object initialized successfully" << std::endl;
//   std::cout << "initial velocities u and v given" << std::endl;
// }
//
// rossby::rossby(double u_in, double v_in, arma::vec psi_init_in)
// {
//   u = u_in;
//   v = v_in;
//   psi_init = psi_init_in;
//   psi.col(0) = psi_init;
//   std::cout << "rossby class object initialized successfully" << std::endl;
//   std::cout << "initial velocities u and v, and initial streamfunction psi given" << std::endl;
// }
//
//  class methods
// void rossby::function()
// {
//   std::cout << "test function for exampleClass, number = " << deltax << std::endl;
// }

//  functions from project 1 for calculating Laplacian

// arma::vec tridiag_general(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int n){
//   // setting up the v vector, giving it the same length as input vector b
//   // b is the diagonal vector
//   arma::vec v(n, arma::fill::zeros);
//   // updating the diagonal vector elements "forward sweep"
//   c(0) = c(0)/b(0);
//   d(0) = d(0)/b(0);
//   for (int i = 1; i <= n-2; i++){
//     c(i) = c(i)/(b(i)-a(i-1)*c(i-1));
//     d(i) = (d(i)-a(i-1)*d(i-1))/(b(i)-a(i-1)*c(i-1));
//   }
//   //backward substitution computes the resulting vector v
//   v(n-1) = (d(n-1)-a(n-2)*d(n-2))/(b(n-1)-a(n-2)*c(n-2)); //v(n-1) = d(n-1)
//   for (int i = n-2; i >= 0; i--) v(i) = d(i)-c(i)*v(i+1);
//   return v;
// }

// Class functions

void rossby::initialize_wave(bool sineWave, double sigma, double x0, double y0)
{
  double x;
  double y;
  for (int i = 0; i < xdim; i++){
    x = (i+1)*deltax;
    for (int j = 0; j < ydim; j++){
      y = (j+1)*deltay;
      if (sineWave)
      {
        double pi = 2*acos(0.0);
        Psi.slice(0)(i,j) = sin(4*pi*x)*sin(4*pi*y);
        Zeta.slice(0)(i,j) = -32*pi*pi*Psi.slice(0)(i,j);
      }
      else
      {
        Psi.slice(0)(i,j) = exp(-pow((x-x0)/sigma,2)-pow((y-y0)/sigma,2));
        Zeta.slice(0)(i,j) = -4/pow(sigma,4)*Psi.slice(0)(i,j)*(sigma*sigma-(x-x0)-(y-y0));
      }
    }
  }
  return;
}

void rossby::zeta_timestep_forward(double &zeta_forward, double zeta, double psi_forward,
  double psi_backward)
{
  zeta_forward = zeta + deltat/(2.0*deltax)*(psi_forward - psi_backward);
  return;
}

void rossby::zeta_timestep_centered(double &zeta_forward, double zeta_backward,
  double psi_forward, double psi_backward)
{
  zeta_forward = zeta_backward + deltat/deltax*(psi_forward - psi_backward);
  return;
}

void rossby::jacobis_method_2d(int n, mat zeta){
  double psiClosed = 0.0;
  double dxdy = deltax*deltay;
  mat psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-6;

  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = Psi.slice(n); difference = 0.;

    //grensenbetingelser sider
    for(int l = 1; l < xdim-1; l++){
      Psi.slice(n)(l,0) = 0.25*(psi_temporary(l,1)+psiClosed
                 +psi_temporary(l+1,0)+psi_temporary(l-1,0)
                 -dxdy*zeta(l,0));
      difference += fabs(psi_temporary(l,0)-Psi.slice(n)(l,0));
      Psi.slice(n)(l,ydim-1) = 0.25*(psiClosed + psi_temporary(l, ydim-2)
                 +psi_temporary(l+1,ydim-1)+psi_temporary(l-1,ydim-1)
                 -dxdy*zeta(l,ydim-1));
      difference += fabs(psi_temporary(l,ydim-1)-Psi.slice(n)(l,ydim-1));
    }

    for(int l = 1; l < ydim-1; l++){
      Psi.slice(n)(0,l) = 0.25*(psi_temporary(0,l+1)+psi_temporary(0,l-1)
                 +psi_temporary(1,l)+psiClosed
                 -dxdy*zeta(0,l));
      difference += fabs(psi_temporary(0,l)-Psi.slice(n)(0,l));
      Psi.slice(n)(xdim-1,l) = 0.25*(psi_temporary(xdim-1,l+1)+psi_temporary(xdim-1,l-1)
                 +psiClosed+psi_temporary(xdim-2,l)
                 -dxdy*zeta(xdim-1,l));
      difference += fabs(psi_temporary(xdim-1,l)-Psi.slice(n)(xdim-1,l));
    }

    //grensebetingelser hjørner
    Psi.slice(n)(0,0) = 0.25*(psi_temporary(0,1)+psiClosed
               +psi_temporary(1,0)+psiClosed
               -dxdy*zeta(0,0));
    difference += fabs(psi_temporary(0,0)-Psi.slice(n)(0,0));
    Psi.slice(n)(xdim-1,ydim-1) = 0.25*(psiClosed+psi_temporary(xdim-1,ydim-2)
               +psiClosed+psi_temporary(xdim-2,ydim-1)
               -dxdy*zeta(xdim-1,ydim-1));
    difference += fabs(psi_temporary(xdim-1,ydim-1)-Psi.slice(n)(xdim-1,ydim-1));
    Psi.slice(n)(0,ydim-1) = 0.25*(psiClosed+psi_temporary(0,ydim-2)
               +psi_temporary(1,ydim-1)+psiClosed
               -dxdy*zeta(0,ydim-1));
    difference += fabs(psi_temporary(0,ydim-1)-Psi.slice(n)(0,ydim-1));
    Psi.slice(n)(xdim-1,0) = 0.25*(psi_temporary(xdim-1,1)+psiClosed
               +psiClosed+psi_temporary(xdim-2,0)
               -dxdy*zeta(xdim-1,0));
    difference += fabs(psi_temporary(xdim-1,0)-Psi.slice(n)(xdim-1,0));

    //ittererer over de indre punktene
    for(int i = 1; i < xdim-1; i++){
      for(int j = 1; j < ydim-1; j++){
        Psi.slice(n)(i,j) = 0.25*(psi_temporary(i,j+1)+psi_temporary(i,j-1)
                   +psi_temporary(i+1,j)+psi_temporary(i-1,j)
                   -dxdy*zeta(i,j));
        difference += fabs(psi_temporary(i,j)-Psi.slice(n)(i,j));
      }
    }
    iterations++;
    difference /= (xdim*ydim);
  }
  return;
}

void rossby::evolve_bounded(bool forwardStep)
{
  double psiClosed = 0.0;
  mat zeta_2previous = Zeta.slice(0);
  mat zeta_previous = Zeta.slice(0);
  for(int n = 0; n < tdim-1; n++){
    for (int j = 0; j < ydim; j++){
      if(forwardStep){
        zeta_timestep_forward(Zeta.slice(n+1)(0,j), zeta_previous(0,j), Psi.slice(n)(1,j), psiClosed);
      }
      else{
        zeta_timestep_centered(Zeta.slice(n+1)(0,j), zeta_2previous(0,j), Psi.slice(n)(1,j), psiClosed);
      }
      for(int i = 1; i < xdim-1; i++){
        if(forwardStep){
          zeta_timestep_forward(Zeta.slice(n+1)(i,j), zeta_previous(i,j), Psi.slice(n)(i+1,j), Psi.slice(n)(i-1,j));
        }
        else{
          zeta_timestep_centered(Zeta.slice(n+1)(i,j), zeta_2previous(i,j), Psi.slice(n)(i+1,j), Psi.slice(n)(i-1,j));
        }
      }
      if(forwardStep){
        zeta_timestep_forward(Zeta.slice(n+1)(xdim-1,j), zeta_previous(xdim-1,j), psiClosed, Psi.slice(n)(xdim-2,j));
      }
      else{
        zeta_timestep_centered(Zeta.slice(n+1)(xdim-1,j), zeta_2previous(xdim-1,j), psiClosed, Psi.slice(n)(xdim-2,j));
      }
    }
    zeta_2previous = zeta_previous;
    zeta_previous = Zeta.slice(n+1);

    jacobis_method_2d(n+1, Zeta.slice(n+1));
  }
  return;
}

void rossby::evolve_periodic(bool forwardStep)
{
  double psiClosed = 0.0;
  mat zeta_2previous = Zeta.slice(0);
  mat zeta_previous = Zeta.slice(0);
  for(int n = 0; n < tdim-1; n++){
    for (int i = 0; i < xdim; i++){
      //Zeta.slice(n+1)(i,0) = psiClosed;
      //for(int j = 1; j < ydim-1; j++){
      for(int j = 0; j < ydim; j++){
        if(forwardStep){
          zeta_timestep_forward(Zeta.slice(n+1)(i,j), zeta_previous(i,j),
          Psi.slice(n)(periodic(i, xdim,1),j), Psi.slice(n)(periodic(i, xdim,-1),j));
        }
        else{
          zeta_timestep_centered(Zeta.slice(n+1)(i,j), zeta_2previous(i,j),
          Psi.slice(n)(periodic(i, xdim,1),j), Psi.slice(n)(periodic(i, xdim,-1),j));
        }
      }
      //Zeta.slice(n+1)(i,ydim-1) = psiClosed;
    }
    zeta_2previous = zeta_previous;
    zeta_previous = Zeta.slice(n+1);

    jacobis_method_2d(n+1, Zeta.slice(n+1));

  }
  return;
}
//
