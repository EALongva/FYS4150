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

Ising::Ising(int N_in, int seed_in)
{
    N = N_in;
    Nspins = N_in*N_in;
    cout << "init dimension = " << N << endl;
    E = 0;
    M = 0;
    J = 1;

    // setting up the e^(-beta*dE) array, hence we dont compute this for
    // every step
    w = arma::vec(17, fill::zeros);
    for (int i=0; i < 17; i += 4){
      w(i) = - (1/T) * (i-8);
    }
    w = exp(w);

    seed = seed_in;

    // setting the arma::random seed
    arma::arma_rng::set_seed(seed);

    // initializing global RNG, define uniform real distribution x \in [[0, 1]
    std::mt19937_64 gen(seed);
    std::uniform_real_distribution<double> URD(0.0,1.0);

    // initialize the state using the Ising::init() function
    init();

}

Ising::Ising(int N_in, double T_in, double J_in, int seed_in, int ordered_in)
{
    N = N_in;
    Nspins = N_in*N_in;
    //cout << "init dimension = " << N << endl;
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

    //cout << "Initial state lattice of dimension = " << N << endl;
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

  // saving the init system
  init_state = state;
  init_M = M;
  init_E = E;

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


void Ising::burninMC(int cycles_in)
{

  int cycles = std::pow(10, cycles_in);

  arma::vec Energies(cycles, arma::fill::zeros);
  arma::vec Magnetisations(cycles, arma::fill::zeros);
  arma::vec AverageEnergies(cycles, arma::fill::zeros);
  arma::vec AverageMagnetisations(cycles, arma::fill::zeros);

  double SumAllEnergies = E/Nspins;
  double SumAllMagnetisations = M/Nspins;

  //std::cout << SumAllEnergies << ", " << SumAllMagnetisations << std::endl;

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

    //std::cout << "sum: " << SumAllEnergies << ", " << SumAllMagnetisations << std::endl;
    //std::cout << "avg: " << AverageEnergies(cyc) << ", " << AverageMagnetisations(cyc) << std::endl;

  }

  avgE = AverageEnergies;
  avgM = AverageMagnetisations;

}


void Ising::MC(int cycles_in)
{

  double cycpow = std::pow(10, cycles_in);
  int cycles = (int) cycpow + 0.5; // + 0.5 to get the correct truncation

  expvals = arma::vec(4, arma::fill::zeros);

  //std::cout << SumAllEnergies << ", " << SumAllMagnetisations << std::endl;

  for (int cyc = 0; cyc < cycles; cyc++){

    for (int spin = 0; spin < Nspins; spin++){

      Metropolis();

    }

    expvals[0] += E;
    expvals[1] += E*E;
    expvals[2] += std::fabs(M);
    expvals[3] += M*M;

  }

  double norm = 1.0/((double) cycles);
  expvals *= norm;

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

  for (int cyc = 0; cyc < cycles; cyc++){

    for (int spin = 0; spin < Nspins; spin++){

      Metropolis();

    }

  }

  burnin_state = state;
  burnin_E = E;
  burnin_M = M;

}


void Ising::observables(int cycles_in, int burnin_cycles_in, double T_start, double T_end, int T_N)
{

  arma::vec T = arma::linspace(T_start, T_end, T_N);

  // values to save

  //double E, E2, M, M2, Mabs;

  for (int k = 0; k < T_N; k++){

    burnin(burnin_cycles_in);
    MC(cycles_in);

  }

  /*
  double norm = 1.0/((double) (MonteCarloCycles));  // divided by  number of cycles
  double E_ExpectationValues = ExpectationValues(0)*norm;
  double E2_ExpectationValues = ExpectationValues(1)*norm;
  double M_ExpectationValues = ExpectationValues(2)*norm;
  double M2_ExpectationValues = ExpectationValues(3)*norm;
  double Mabs_ExpectationValues = ExpectationValues(4)*norm; */

}

void Ising::reset2burnin(){
  if (burnin_state[0,0] == 0){
    std::cout << "no burn in state exist" << std::endl;
  }
  else{
    state = burnin_state;
    E = burnin_E;
    M = burnin_M;
  }
}

void Ising::reset2init(){

    state = init_state;
    E = init_E;
    M = init_M;

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


void Ising::hello(int &argc, char** &argv){

  int NP, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::cout << "Rank: " << rank << ", NP: " << NP << std::endl;

  MPI_Finalize();

}

void Ising::parallel(int &argc, char** &argv){

  // declaring variables
  int burnin_cycles, MC_cycles, Tlen;
  double T_init, T_fin, dT;
  std::string filename;
  std::string outpath = "results/data/";

  int NP, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &NP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) {
    if (argc > 1){
      filename   = outpath+argv[1];
      N          = atoi(argv[2]);
      burnin_cycles     = atoi(argv[3]);
      MC_cycles         = atoi(argv[4]);
      T_init     = atof(argv[5]);
      T_fin      = atof(argv[6]);
      dT         = atof(argv[7]);

      Tlen       = (int) ((T_fin - T_init)/dT + 0.5);
      totexpvals    = arma::mat(4,Tlen,arma::fill::zeros);

    }
    else{
      // standard values (aka lazy values)
      N               = 20;
      burnin_cycles   = 5;
      MC_cycles       = 5;

      T_init     = 2.0;
      T_fin      = 2.4;
      dT         = 0.05;

      Tlen        = (int) ((T_fin - T_init)/dT + 0.5) + 1;
      totexpvals  = arma::mat(4,Tlen,arma::fill::zeros);

      filename = outpath + "Simulation_N" + std::to_string(N) + ".csv";

      std::cout << "No cmd line args given: " << argv[0] << std::endl;
      std::cout << " read output file, Number of spins, burn-in MC cycles, MC cycles,"
      << "initial and final temperature and tempurate step" << std::endl;
      std::cout << "-------------------------------------------------" << std::endl;
      std::cout << "Implementing standard values: " << filename << ", " << N <<
      ", " << burnin_cycles << ", " << MC_cycles << ", " << T_init << ", " << T_fin << ", "<<
      dT << std::endl;
    }
  }

  // broadcast to all nodes common variables since only master node reads from command line
  MPI_Bcast (&burnin_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&MC_cycles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&T_init, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&T_fin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // fix above variables and start the main temperature loop while timing it
  // check the amount of cycles run

  double t_init, t_fin, t_total;
  t_init = MPI_Wtime();

  int i = 0;

  for (double temp = T_init; temp <= T_fin; temp += dT){

    T = temp;
    init();
    burnin(burnin_cycles);

    MC(MC_cycles);

    for (int j=0; j<4; j++){
      //std::cout << "j: " << j << ", i: " << i << ", Tlen:" << Tlen << std::endl;
      //MPI_Reduce(&expvals(j), &totexpvals(j,i), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    i ++;
  }

  if (rank == 0){
    totexpvals *= 1/((double) NP);
    totexpvals.save(filename, arma::csv_ascii);
  }

  t_fin = MPI_Wtime();
  t_total = t_fin - t_init;

  if (rank == 0){
    std::cout << "total time: " << t_total << ", using: " << NP << "processors" << std::endl;
    std::cout << "number of Monte Carlo cycles = " << std::pow(10,MC_cycles)*NP << std::endl;
  }

  MPI_Finalize();

}

//
