#include "../include/metro.hpp"
using namespace arma;
using namespace std;

// inline int periodic(int i, int limit, int add){
//   return (i+limit+add) % (limit);}

void set_spin(double& val)
{
  if (val <= 0.5){ val = -1; }
  else { val = 1; }
}

Ising::Ising(int N_in)
{
    N = N_in;
    cout << "init dimension = " << N << endl;
}

Ising::Ising(int N_in, double T_in, double J_in, int seed)
{
    N = N_in;
    cout << "init dimension = " << N << endl;

    T = T_in;
    J = J_in;
    arma_rng::set_seed(seed);

}

void Ising::genState()
{
  arma_rng::set_seed_random();
  state.randu(N,N);
  state.for_each( [](arma::mat::elem_type& val) { set_spin(val); } );
  state.col(0).print();

  cout << "Initial state lattice of dimension = " << N << endl;
  state.print();
}


void Ising::energy()
{
  double tempE = 0;
  mat pS; // altered state having periodic boundaries (periodic State)

  // concatenation of upper and lower rows to be looped over
  // -> a somewhat tricky process, appends the last column to the
  //   start of the matrix, and the first original column to the
  //   end of the matrix, same procedure for rows
  //   this gives a matrix with periodic boundaries

  pS = join_rows(state.col(N-1), state);
  pS = join_rows(pS, state.col(0));
  pS = join_cols(pS.row(N-1), pS);
  pS = join_cols(pS, pS.row(0));

  // pS.print();
  // std::cout << pS(0,0) << std::endl;

  for (int i = 1; i <= N; i++){

    for (int j = 1; j <= N; j++){

      tempE += -J*( pS(i,j)*pS(i-1,j) + pS(i,j)*pS(i+1,j)
                    + pS(i,j)*pS(i,j-1) + pS(i,j)*pS(i,j+1) );

    }

    cout << tempE << endl;

  }

}

void Ising::Metropolis()
{
  for (int y = 0; y < L; y++){
    for (int x = 0; y < L; x++){
      int ix = (int) rand() % L;
      int iy = (int) rand() % L;
      int deltaE = 2 * spin_state(iy, ix) *
    }
  }

  // void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w)
  // {
  //   // loop over all spins
  //   for(int y =0; y < n_spins; y++) {
  //     for (int x= 0; x < n_spins; x++){
  //       int ix = (int) (ran1(&idum)*(double)n_spins);
  //       int iy = (int) (ran1(&idum)*(double)n_spins);
  //       int deltaE =  2*spin_matrix[iy][ix]*
  // 	(spin_matrix[iy][periodic(ix,n_spins,-1)]+
  // 	 spin_matrix[periodic(iy,n_spins,-1)][ix] +
  // 	 spin_matrix[iy][periodic(ix,n_spins,1)] +
  // 	 spin_matrix[periodic(iy,n_spins,1)][ix]);
  //       if ( ran1(&idum) <= w[deltaE+8] ) {
  // 	spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
  //         M += (double) 2*spin_matrix[iy][ix];
  //         E += (double) deltaE;
  //       }
  //     }
  //   }
  // } // end of Metropolis sampling over spins


}


//
