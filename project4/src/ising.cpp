#include "../include/ising.hpp"
using namespace arma;
using namespace std;

inline int periodic(int i, int limit, int add){
  return (i+limit+add) % (limit);}

void set_spin(double& val)
{
  if (val <= 0.5){ val = ((int) -1); }
  else { val = ((int) 1); }
}

Ising::Ising(int N_in, double T_in)
{
    N = N_in;
    cout << "init dimension = " << N << endl;
    E = 0;
    M = 0;
    T = T_in;
    J = 1;

    // setting up the e^(-beta*dE) array, hence we dont compute this for
    // every step
    w = arma::vec(17, fill::zeros);
    for (int i=0; i < 17; i += 4){
      w(i) = - (1/T) * (i-8);
    }
    w = exp(w);
}

Ising::Ising(int N_in, double T_in, double J_in, int seed_in, int ordered_in)
{
    N = N_in;
    Nspins = N_in*N_in;
    cout << "init dimension = " << N << endl;
    E = 0;
    M = 0;
    T = T_in;
    J = J_in;

    // setting up the e^(-beta*dE) array, hence we dont compute this for
    // every step
    w = arma::vec(17, fill::zeros);
    for (int i=0; i < 17; i += 4){
      w(i) = - (1/T) * (i-8);
    }
    w = exp(w);

    seed = seed_in;
    ordered = ordered_in;

    // setting the arma::random seed
    arma::arma_rng::set_seed(seed);

    // initializing global RNG, define uniform real distribution x \in [[0, 1]
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> URD(0.0,1.0);

    // initialize the state using the Ising::init() function
    init();

}

void Ising::init()
{
  if (T <= 2.1 or ordered == 0){

    state = arma::mat(N, N, arma::fill::ones);

  }

  else{

    arma_rng::set_seed(seed);

    state.randu(N,N);
    state.for_each( [](arma::mat::elem_type& val) { set_spin(val); } );

    cout << "Initial state lattice of dimension = " << N << endl;
    state.print();
  }

  double tempE = 0;
  mat pS; // altered state having periodic boundaries (periodic State)

  // description in ising::energy()

  pS = join_rows(state.col(N-1), state);
  pS = join_rows(pS, state.col(0));
  pS = join_cols(pS.row(N-1), pS);
  pS = join_cols(pS, pS.row(0));

  for (int i = 1; i <= N; i++){
    for (int j = 1; j <= N; j++){

      tempE += -J*pS(i,j)*( pS(i-1,j) + pS(i,j-1) );

    }
  }

  E = tempE;
  M = arma::accu(state);

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
  // -> a some what inefficient process, appends the last column to the
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
  // Performs one metropolis alogorithm

  int ix = (int) (URD(gen)*N);
  int iy = (int) (URD(gen)*N);

  /*std::cout << "random indices: " << ix << iy << std::endl;
  std::cout << "periodic neighbors: " << std::endl; */

  // Computing the energy change dE
  int dE = 2 * state(iy, ix) *
    (state(iy, periodic(ix, N, -1)) +
    state(periodic(iy, N, -1), ix) +
    state(iy, periodic(ix, N, 1)) +
    state(periodic(iy, N, 1), ix));

  // We use our energy change dE to index a corresponding probability w,
  // which is dependent on temperature. A dE=-8 corresponds to
  // a probability w(dE+8)=w(0)=exp(8/(k_bT)).
  // We flip the spin if the energy change is less than or equal to
  // zero (which gives a w(dE+8)=>1) or if a generated random number between
  // 0 and 1 is less than w(dE+8)

  if (URD(gen) <= w(dE+8)){
    state(iy, ix) *= -1; // flip the spin
    M += 2*state(iy, ix);
    E += dE;
  }
}

void Ising::MC(int cycles_in, int burnin)
{
  if(burnin > 1 or burnin < 0){
    std::cout << "invalid value for burnin, 0 = burnin false, 1 = burnin true" << std::endl;
  }

  int cycles = std::pow(10, cycles_in);

  arma::vec Energies(cycles, arma::fill::zeros);
  arma::vec Magnetisations(cycles, arma::fill::zeros);
  arma::vec AverageEnergies(cycles, arma::fill::zeros);
  arma::vec AverageMagnetisations(cycles, arma::fill::zeros);

  double SumAllEnergies = E/Nspins;
  double SumAllMagnetisations = M/Nspins;

  std::cout << SumAllEnergies << ", " << SumAllMagnetisations << std::endl;

  for (int cyc = 0; cyc < cycles; cyc++){

    for (int spin = 0; spin < Nspins; spin++){

      Metropolis();

    }

    Energies(cyc) = E/Nspins;
    Magnetisations(cyc) = fabs(M)/Nspins;
    SumAllEnergies += Energies(cyc);
    SumAllMagnetisations += Magnetisations(cyc);
    AverageEnergies(cyc) = SumAllEnergies/(cyc+1);
    AverageMagnetisations(cyc) = SumAllMagnetisations/(cyc+1);

    std::cout << "sum: " << SumAllEnergies << ", " << SumAllMagnetisations << std::endl;
    std::cout << "avg: " << AverageEnergies(cyc) << ", " << AverageMagnetisations(cyc) << std::endl;

  }

  avgE = AverageEnergies;
  avgM = AverageMagnetisations;

}


void Ising::paraMC(int cycles_in)
{

  // make function that runs standard Monte Carlo in parallel, loop over time
  // in other function

  int cycles = std::pow(10, cycles_in);

  arma::vec Energies(cycles, arma::fill::zeros);
  arma::vec Magnetisations(cycles, arma::fill::zeros);
  arma::vec AverageEnergies(cycles, arma::fill::zeros);
  arma::vec AverageMagnetisations(cycles, arma::fill::zeros);

  double SumAllEnergies = E/Nspins;
  double SumAllMagnetisations = M/Nspins;

  for (int cyc = 0; cyc < cycles; cyc++){

    for (int spin = 0; spin < Nspins; spin++){

      Metropolis();

    }

    Energies(cyc) = E/Nspins;
    Magnetisations(cyc) = fabs(M)/Nspins;
    SumAllEnergies += Energies(cyc);
    SumAllMagnetisations += Magnetisations(cyc);
    AverageEnergies(cyc) = SumAllEnergies/(cyc+1);
    AverageMagnetisations(cyc) = SumAllMagnetisations/(cyc+1);

  }

  avgE = AverageEnergies;
  avgM = AverageMagnetisations;

}


void Ising::burnin(int cycles_in)
{

  // cycles_in how many cycles to burn

  int cycles = std::pow(10, cycles_in);

  arma::vec Energies(cycles, arma::fill::zeros);
  arma::vec Magnetisations(cycles, arma::fill::zeros);
  arma::vec AverageEnergies(cycles, arma::fill::zeros);
  arma::vec AverageMagnetisations(cycles, arma::fill::zeros);

  double SumAllEnergies = E/Nspins;
  double SumAllMagnetisations = M/Nspins;

  for (int cyc = 0; cyc < cycles; cyc++){

    for (int spin = 0; spin < Nspins; spin++){

      Metropolis();

    }

    Energies(cyc) = E/Nspins;
    Magnetisations(cyc) = fabs(M)/Nspins;
    SumAllEnergies += Energies(cyc);
    SumAllMagnetisations += Magnetisations(cyc);
    AverageEnergies(cyc) = SumAllEnergies/(cyc+1);
    AverageMagnetisations(cyc) = SumAllMagnetisations/(cyc+1);

  }

  avgE = AverageEnergies;
  avgM = AverageMagnetisations;

}


void Ising::save(std::string filename)
{
  // saves average energies and magnetisations using the arma::save()

  int length          = avgE.n_elem;
  std::string len     = std::to_string(length);
  std::string csv     = ".csv";

  std::string fileout_avgE = "results/data/" + filename + "_avgE_" + len + csv;
  std::string fileout_avgM = "results/data/" + filename + "_avgM_" + len + csv;

  avgE.save(fileout_avgE, arma::csv_ascii);
  avgM.save(fileout_avgM, arma::csv_ascii);

}


//
