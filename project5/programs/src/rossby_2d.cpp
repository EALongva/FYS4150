#include "../include/rossby_2d.hpp"

using namespace arma;
using namespace std;

// Modulo indexing to get periodic boundary conditions
inline int periodic(int i, int limit, int add){
  return (i+limit+add) % (limit);}


//  Initializer

rossby::rossby(double dpos, double dt, double tfinal)
{
  endposx = 1.0; // Length of spatial domain in x-direction
  endposy = 1.0; // Length of spatial domain in y-direction
  endtime = tfinal; // Length of time period
  xdim = (int) endposx/dpos; // Spatial x dimension
  ydim = (int) endposy/dpos; // Spatial y dimension
  tdim = (int) endtime/dt; // Temporal dimension
  deltax = dpos; // x grid space
  deltay = dpos; // y grid space
  deltat = dt; // time step
  Psi = cube(xdim, ydim, tdim, fill::zeros); // Streamfunction matrix
  Zeta = cube(xdim, ydim, tdim, fill::zeros); // Vorticity cube
  cout << "rossby class object initialized successfully" << endl;
}

// Class functions

void rossby::initialize_wave(bool sineWave, double sigma, double x0, double y0)
{
  // Initialize wave in first time point, either a sine or a Gaussian
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
        Zeta.slice(0)(i,j) = -4/pow(sigma,4)*Psi.slice(0)(i,j)*(sigma*sigma-pow(x-x0,2)-pow(y-y0,2));
      }
    }
  }
  return;
}

void rossby::zeta_timestep_forward(double &zeta_forward, double zeta, double psi_forward,
  double psi_backward)
{
  // Forward difference time step
  zeta_forward = zeta - deltat/(2.0*deltax)*(psi_forward - psi_backward);
  return;
}

void rossby::zeta_timestep_centered(double &zeta_forward, double zeta_backward,
  double psi_forward, double psi_backward)
{
  // Centered difference timestep
  zeta_forward = zeta_backward - deltat/deltax*(psi_forward - psi_backward);
  return;
}

void rossby::jacobis_method_2d_bounded(int n, mat zeta){

  double psiClosed = 0.0;
  double dxdy = deltax*deltay;
  mat psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-6;

  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = Psi.slice(n); difference = 0.;

    // Interate over all edges
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

    // Iterate over each corner
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

    // Iterate over interior points
    for(int i = 1; i < xdim-1; i++){
      for(int j = 1; j < ydim-1; j++){
        Psi.slice(n)(i,j) = 0.25*(psi_temporary(i,j+1)+psi_temporary(i,j-1)
                   +psi_temporary(i+1,j)+psi_temporary(i-1,j)
                   -dxdy*zeta(i,j));
        difference += fabs(psi_temporary(i,j)-Psi.slice(n)(i,j));
      }
    }
    iterations++;
    difference /= (xdim*ydim); // Divide difference by number of points
  }
  return;
}

void rossby::jacobis_method_2d_periodic(int n, mat zeta){

  double psiClosed = 0.0;
  double dxdy = deltax*deltay;
  mat psi_temporary;
  int iterations = 0; int maxIterations = 10000;
  double difference = 1.; double maxDifference = 1e-6;

  while((iterations <= maxIterations) && (difference > maxDifference)){
    psi_temporary = Psi.slice(n); difference = 0.;

    // Interate over along northern and southern edges
    for(int l = 0; l < xdim; l++){
      Psi.slice(n)(l,0) = 0.25*(psi_temporary(l,1)+psiClosed
                 +psi_temporary(periodic(l, xdim,1),0)+psi_temporary(periodic(l, xdim,-1),0)
                 -dxdy*zeta(l,0));
      difference += fabs(psi_temporary(l,0)-Psi.slice(n)(l,0));
      Psi.slice(n)(l,ydim-1) = 0.25*(psiClosed + psi_temporary(l, ydim-2)
                 +psi_temporary(periodic(l, xdim,1),ydim-1)+psi_temporary(periodic(l, xdim,-1),ydim-1)
                 -dxdy*zeta(l,ydim-1));
      difference += fabs(psi_temporary(l,ydim-1)-Psi.slice(n)(l,ydim-1));
    }

    // Iterate over interior points
    for(int i = 0; i < xdim; i++){
      for(int j = 1; j < ydim-1; j++){
        Psi.slice(n)(i,j) = 0.25*(psi_temporary(i,j+1)+psi_temporary(i,j-1)
                   +psi_temporary(periodic(i, xdim,1),j)+psi_temporary(periodic(i, xdim,-1),j)
                   -dxdy*zeta(i,j));
        difference += fabs(psi_temporary(i,j)-Psi.slice(n)(i,j));
      }
    }
    iterations++;
    difference /= (xdim*ydim); // Divide difference by number of points
  }
  return;
}


void rossby::evolve_bounded(bool forwardStep)
{
  double psiClosed = 0.0;
  mat zeta_2previous = Zeta.slice(0);
  mat zeta_previous = Zeta.slice(0);
  // looping over each time step
  for(int n = 0; n < tdim-1; n++){
    // calculating the new voriticty at the left boundary
    for (int j = 0; j < ydim; j++){
      if(forwardStep){
        zeta_timestep_forward(Zeta.slice(n+1)(0,j), zeta_previous(0,j), Psi.slice(n)(1,j), psiClosed);
      }
      else{
        zeta_timestep_centered(Zeta.slice(n+1)(0,j), zeta_2previous(0,j), Psi.slice(n)(1,j), psiClosed);
      }
        // calculating the new interior vorticity
      for(int i = 1; i < xdim-1; i++){
        if(forwardStep){
          zeta_timestep_forward(Zeta.slice(n+1)(i,j), zeta_previous(i,j), Psi.slice(n)(i+1,j), Psi.slice(n)(i-1,j));
        }
        else{
          zeta_timestep_centered(Zeta.slice(n+1)(i,j), zeta_2previous(i,j), Psi.slice(n)(i+1,j), Psi.slice(n)(i-1,j));
        }
      }
      // calculating the new vorticity at the right boundary
      if(forwardStep){
        zeta_timestep_forward(Zeta.slice(n+1)(xdim-1,j), zeta_previous(xdim-1,j), psiClosed, Psi.slice(n)(xdim-2,j));
      }
      else{
        zeta_timestep_centered(Zeta.slice(n+1)(xdim-1,j), zeta_2previous(xdim-1,j), psiClosed, Psi.slice(n)(xdim-2,j));
      }
    }
    // storing the two previous vortitices two use in next timestep
    zeta_2previous = zeta_previous;
    zeta_previous = Zeta.slice(n+1);

    // updating the streamfunction
    jacobis_method_2d_bounded(n+1, Zeta.slice(n+1));
  }
  return;
}

void rossby::evolve_periodic(bool forwardStep)
{
  double psiClosed = 0.0;
  mat zeta_2previous = Zeta.slice(0);
  mat zeta_previous = Zeta.slice(0);

  for(int n = 0; n < tdim-1; n++){
    // calculating the new vorticity in each spatial point using periodic boundaries
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
    // storing the two previous vortitices two use in next timestep
    zeta_2previous = zeta_previous;
    zeta_previous = Zeta.slice(n+1);

    // updating the streamfunction
    jacobis_method_2d_periodic(n+1, Zeta.slice(n+1));

  }
  return;
}
//
