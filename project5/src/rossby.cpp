#include "../include/rossby.hpp"

using namespace arma;
using namespace std;

inline int periodic(int i, int limit, int add){
  return (i+limit+add) % (limit);}


//  Initializers

rossby::rossby(double dx, double dt, double tfinal)
{
  endpos = 1.0;
  endtime = tfinal;
  xdim = (int) endpos/dx; 
  tdim = (int) endtime/dt;
  deltax = dx;
  deltat = dt;
  Psi = mat(xdim, tdim, fill::zeros);
  Zeta = mat(xdim, tdim, fill::zeros);
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

void rossby::initialize_wave(bool sineWave, double sigma, double x0)
{
  double x;
  for (int j = 0; j < xdim; j++){
    x = (j+1)*deltax;
    if (sineWave)
    {
      double pi = 2*acos(0.0);
      Psi.col(0)(j) = sin(4*pi*x);
      Zeta.col(0)(j) = -16*pi*pi*sin(4*pi*x);
    }
    else
    {
      Psi.col(0)(j) = exp(-pow((x-x0)/sigma,2));
      Zeta.col(0)(j) = -2/pow(sigma,4)*Psi.col(0)(j)*(sigma*sigma - 2*(x-x0));
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

vec rossby::precalculate_offdiag()
{
  vec c_new(xdim, fill::zeros);
  for (int j = 1; j < xdim+1; j++){
    c_new(j-1) = - (double) j / (j+1);
  }
  return c_new;
}

void rossby::gaussian_elimination(int n, vec c_new)
{
  vec d = Zeta.col(n)*deltax*deltax;
  vec d_new(xdim, fill::zeros);

  // forward substitution to find the new equations
  d_new(0) = d(0)*c_new(0);
  for (int j = 1; j < xdim; j++){
    d_new(j) = c_new(j)*(d(j) - d_new(j-1));
  }

  // backward substitution to find the new psi
  Psi.col(n)(xdim-1) = d_new(xdim-1);
  for (int j = xdim-2; j >= 0; j--){
    Psi.col(n)(j) = d_new(j) - c_new(j)*Psi.col(n)(j+1);
  }
  return;
}

//
// vec rossby::forward_substitution(vec zeta, vec c_new){
//   vec d = zeta*deltax*deltax;
//   vec d_new(xdim, fill::zeros);
//
//   // forward substitution to find the new equation
//   d_new(0) = d(0)*c_new(0);
//   for (int j = 1; j < xdim; j++){
//     d_new(j) = c_new(j)*(d(j) - d_new(j-1));
//   }
//   return d_new;
// }

void rossby::jacobis_method(int n, vec zeta){
  double hh = deltax*deltax;
  vec psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-6;
  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = Psi.col(n); difference = 0.;
    for(int j = 0; j < xdim; j++){
      Psi.col(n)(j) = 0.5*(psi_temporary(periodic(j, xdim,1)) + psi_temporary(periodic(j, xdim,-1)) - zeta(j)*hh);
      difference += fabs(psi_temporary(j)-Psi.col(n)(j));
    }
    iterations++;
    difference /= xdim;
  }
  return;
}
// arma::vec tridiag_LU(arma::mat A, arma::vec d, int n){
//
//   // setting up the v vector, giving it the same length as input vector b
//   // v is the solution vector (u in project description)
//   arma::vec v(n, arma::fill::zeros);
//
//   // also setup the intermidiate step vector y to solve the LU problem
//   arma::vec y(n, arma::fill::zeros);
//
//   // using the LU decomp from armadillo
//   arma::mat L;
//   arma::mat U;
//   arma::mat P;
//   arma::lu(L, U, P, A);
//
//   // using arma::solve to solve the linear problem (2 steps)
//   d = P * d; // unnecessary for our matrix ?
//   y = arma::solve(L,d);
//   v = arma::solve(U,y);
//
//   return v;
// }


void rossby::evolve_bounded(bool forwardStep)
{
  double psiClosed = 0;
  vec zeta_2previous = Zeta.col(0);
  vec zeta_previous = Zeta.col(0);
  vec c_new = precalculate_offdiag();
  for(int n = 0; n < tdim-1; n++){
    if(forwardStep){
      zeta_timestep_forward(Zeta.col(n+1)(0), zeta_previous(0), Psi.col(n)(1), psiClosed);
    }
    else{
      zeta_timestep_centered(Zeta.col(n+1)(0), zeta_2previous(0), Psi.col(n)(1), psiClosed);
    }
    for(int j = 1; j < xdim-1; j++){
      if(forwardStep){
        zeta_timestep_forward(Zeta.col(n+1)(j), zeta_previous(j), Psi.col(n)(j+1), Psi.col(n)(j-1));
      }
      else{
        zeta_timestep_centered(Zeta.col(n+1)(j), zeta_2previous(j), Psi.col(n)(j+1), Psi.col(n)(j-1));
      }
    }
    if(forwardStep){
      zeta_timestep_forward(Zeta.col(n+1)(xdim-1), zeta_previous(xdim-1), psiClosed, Psi.col(n)(xdim-2));
    }
    else{
      zeta_timestep_centered(Zeta.col(n+1)(xdim-1), zeta_2previous(xdim-1), psiClosed, Psi.col(n)(xdim-2));
    }
    zeta_2previous = zeta_previous;
    zeta_previous = Zeta.col(n+1);

    gaussian_elimination(n+1, c_new);

    // // Gaussian eliminiation
    // // Forward substitution
    // vec d_new = forward_substitution(Zeta.col(n+1),c_new);
    //
    // // Backward substitution
    // Psi(xdim-1,n+1) = d_new(xdim-1);
    // for (int j = xdim-2; j = 0; j--){
    //   Psi(j,n+1) = d_new(j) - c_new(j)*Psi(j+1,n+1);
    // }
  }
  return;
}

void rossby::evolve_periodic(bool forwardStep)
{
  vec zeta_2previous = Zeta.col(0);
  vec zeta_previous = Zeta.col(0);
  for(int n = 0; n < tdim-1; n++){
    // finner den første x-verdien til zeta
    for(int j = 0; j < xdim; j++){
      if(forwardStep){
        zeta_timestep_forward(Zeta.col(n+1)(j), zeta_previous(j), Psi.col(n)(periodic(j, xdim,1)), Psi.col(n)(periodic(j, xdim,-1)));
      }
      else{
        zeta_timestep_centered(Zeta.col(n+1)(j), zeta_2previous(j), Psi.col(n)(periodic(j, xdim,1)), Psi.col(n)(periodic(j, xdim,-1)));
      }
    }
    zeta_2previous = zeta_previous;
    zeta_previous = Zeta.col(n+1);

    jacobis_method(n+1, Zeta.col(n+1));

  }
  return;
}
//
