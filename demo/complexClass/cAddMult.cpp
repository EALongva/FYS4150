// complex addition and multiplication

# include <iostream>
# include <cstdlib>
# include <cmath>
# include <complex>

using namespace std;

int main (int argc, char* argv[])
{
  cout << "number of arguments given =   " << argc - 1 << endl;
  // exception if argc =/= 5
  if (argc != 5){
    cout << "incorrect number of arguments" << endl;
    cout << "please insert 4 args to represent two complex numbers" << endl;
    exit(1);
  }

  double a1 = atof(argv[1]);
  double b1 = atof(argv[2]);
  double a2 = atof(argv[3]);
  double b2 = atof(argv[4]);

  complex<double> z1(a1, b1), z2(a2, b2);

  cout << "adding two complex numbers" << z1 + z2 << endl;
  cout << "multiplying two complex numbers" << z1*z2 << endl;

  return 0;
}
