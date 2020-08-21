// creating a vector

# include <iostream>
# include <armadillo>
using namespace std;
using namespace arma;

int main(){
  A = mat(4, 4); // matrix size (rows,cols)

  A << 1 << 2 << 3 << 4 << endr
    << 2 << 1 << 2 << 3 << endr
    << 3 << 2 << 1 << 2 << endr
    << 4 << 3 << 2 << 1 << endr;

  A.print("A:");

  cout << " The rank of A is " << arma::rank(A) << endl;

  return 0;
}
