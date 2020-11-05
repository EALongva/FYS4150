#include "../include/metro.hpp"

//using namespace arma;
//using namespace std;

void set_spin(double& val){

  if (val <= 0.5){ val = -1; }
  else { val = 1; }

}

metro::metro(int N_in)
{
    N = N_in;
    std::cout << "init dimension = " << N << std::endl;
}

metro::metro(int N_in, double T_in, double J_in, int seed)
{
    N = N_in;
    std::cout << "init dimension = " << N << std::endl;

    T = T_in;
    J = J_in;
    arma::arma_rng::set_seed(seed);

}

void metro::genState()
{
  arma::arma_rng::set_seed_random();
  state.randu(N,N);
  state.for_each( [](arma::mat::elem_type& val) { set_spin(val); } );
  state.col(0).print();

  std::cout << "Initial state lattice of dimension = " << N << std::endl;
  state.print();
}


void metro::energy()
{
  double tempE = 0;
  arma::mat pS; // altered state having periodic boundaries (periodic State)

  // concatenation of upper and lower rows to be looped over
  // -> a somewhat tricky process, appends the last column to the
  //   start of the matrix, and the first original column to the
  //   end of the matrix, same procedure for rows
  //   this gives a matrix with periodic boundaries

  pS = arma::join_rows(state.col(N-1), state);
  pS = arma::join_rows(pS, state.col(0));
  pS = arma::join_cols(pS.row(N-1), pS);
  pS = arma::join_cols(pS, pS.row(0));

  // pS.print();
  // std::cout << pS(0,0) << std::endl;

  for (int i = 1; i <= N; i++){

    for (int j = 1; j <= N; j++){

      tempE += -J*( pS(i,j)*pS(i-1,j) + pS(i,j)*pS(i+1,j)
                    + pS(i,j)*pS(i,j-1) + pS(i,j)*pS(i,j+1) );

    }

    std::cout << tempE << std::endl;

  }

}

//
