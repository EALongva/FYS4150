#include "../include/rossby.hpp"

using namespace arma;
using namespace std;

//  Initializers

rossby::rossby(arma::vec psi_init_in)
{
  psi_init = psi_init_in;
  psi.col(0) = psi_init;
  std::cout << "rossby class object initialized successfully" << std::endl;
  std::cout << "initial streamfunction psi given" << std::endl;
}

rossby::rossby(double u_in, double v_in)
{
  u = u_in;
  v = v_in;
  std::cout << "rossby class object initialized successfully" << std::endl;
  std::cout << "initial velocities u and v given" << std::endl;
}

rossby::rossby(double u_in, double v_in, arma::vec psi_init_in)
{
  u = u_in;
  v = v_in;
  psi_init = psi_init_in;
  psi.col(0) = psi_init;
  std::cout << "rossby class object initialized successfully" << std::endl;
  std::cout << "initial velocities u and v, and initial streamfunction psi given" << std::endl;
}

//  class methods
void rossby::function()
{
  std::cout << "test function for exampleClass, number = " << number << std::endl;
}

//  functions from project 1 for calculating Laplacian

arma::vec tridiag_general(arma::vec a, arma::vec b, arma::vec c, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n, arma::fill::zeros);
  // updating the diagonal vector elements "forward sweep"
  c(0) = c(0)/b(0);
  d(0) = d(0)/b(0);
  for (int i = 1; i <= n-2; i++){
    c(i) = c(i)/(b(i)-a(i-1)*c(i-1));
    d(i) = (d(i)-a(i-1)*d(i-1))/(b(i)-a(i-1)*c(i-1));
  }
  //backward substitution computes the resulting vector v
  v(n-1) = (d(n-1)-a(n-2)*d(n-2))/(b(n-1)-a(n-2)*c(n-2)); //v(n-1) = d(n-1)
  for (int i = n-2; i >= 0; i--) v(i) = d(i)-c(i)*v(i+1);
  return v;
}

arma::vec gaussian_eliminiation(arma::vec b, arma::vec d, int n){
  // setting up the v vector, giving it the same length as input vector b
  // b is the diagonal vector
  arma::vec v(n, arma::fill::zeros);
  // updating the diagonal vector elements "forward sweep" (not looping over d(n-1))
  for (int i = 1; i < n+1; i++){
    b(i-1) = (i+1.0)/((double) i);
  }
  for (int i = 1; i < n; i++){
    d(i) = d(i) + d(i-1)/b(i-1);
  }
  //backward substitution computes the resulting vector v
  v(n-1) = d(n-1)/b(n-1);
  for (int i = n-2; i >= 0; i--) v(i) = (d(i) + v(i+1))/b(i);
  return v;
}

arma::vec tridiag_LU(arma::mat A, arma::vec d, int n){

  // setting up the v vector, giving it the same length as input vector b
  // v is the solution vector (u in project description)
  arma::vec v(n, arma::fill::zeros);

  // also setup the intermidiate step vector y to solve the LU problem
  arma::vec y(n, arma::fill::zeros);

  // using the LU decomp from armadillo
  arma::mat L;
  arma::mat U;
  arma::mat P;
  arma::lu(L, U, P, A);

  // using arma::solve to solve the linear problem (2 steps)
  d = P * d; // unnecessary for our matrix ?
  y = arma::solve(L,d);
  v = arma::solve(U,y);

  return v;
}

void timestep_forward(double &zeta_forward, double psi_forward, double psi_backward, double deltat, double deltax)
{
  zeta_forward += deltat/(2*deltax)*(psi_forward - psi_backward)
}

void timestep_centered(double &zeta_forward, double zeta_backward, double psi_forward, double psi_backward, double deltat, double deltax)
{
  zeta_forward = zeta_backward + deltat/deltax*(psi_forward - psi_backward))
}



//
