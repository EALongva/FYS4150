#include "../include/rossby.hpp"

using namespace arma;
using namespace std;

//  Initializers

rossby::rossby(double dx, double dt, int J, int N)
{
  xdim = J;
  tdim = N;
  deltax = dx;
  deltat = dt;
  psi = mat(xdim, tdim, fill::zeros);
  zeta = mat(xdim, tdim, fill::zeros);
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
void rossby::function()
{
  std::cout << "test function for exampleClass, number = " << deltax << std::endl;
}

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

void rossby::initialize_wave(bool sineWave, double sigma, double x0)
{
  double x;
  for (int j = 0; j < xdim; j++){
    x = (j+1)*deltax;
    if (sineWave)
    {
      double pi = 2*acos(0.0);
      psi.col(0)(j) = sin(4*pi*x);
      zeta.col(0)(j) = -16*pi*pi*sin(4*pi*x);
    }
    else
    {
      psi.col(0)(j) = exp(-pow((x-x0)/sigma,2));
      zeta.col(0)(j) = -2/pow(sigma,4)*psi(j)*(sigma*sigma - 2*(x-x0));
    }
  }
  return;
}


void rossby::precalculate_offdiag(vec &c_new)
{
  for (int j = 1; j < xdim+1; j++){
    c_new(j-1) = - (double) j / (j+1);
  }
  return;
}

void rossby::gaussian_elimination(vec &psi_forward, vec zeta, vec c_new){
  vec d = zeta*deltax*deltax;
  vec d_new(xdim, fill::zeros);

  // forward substitution to find the new equation
  d_new(0) = d(0)*c_new(0);
  for (int j = 1; j < xdim; j++){
    d_new(j) = c_new(j)*(d(j) - d_new(j-1));
  }

  // backward substitution computes the resulting vector psi_forward
  psi_forward(xdim-1) = d_new(xdim-1);
  for (int j = xdim-2; j = 0; j--){
    psi_forward(j) = d_new(j) - c_new(j)*psi_forward(j+1);
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

void rossby::timestep_forward(double &zeta_forward, double psi_forward,
  double psi_backward)
{
  zeta_forward += deltat/(2*deltax)*(psi_forward - psi_backward);
  return;
}

void rossby::timestep_centered(double &zeta_forward, double zeta_backward,
  double psi_forward, double psi_backward)
{
  zeta_forward = zeta_backward + deltat/deltax*(psi_forward - psi_backward);
  return;
}



//
