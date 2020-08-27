// cmain, making complex computations with the Complex Class

#include <Complex.h>
#include <iostream>
using namespace std;

int main()
{
  Complex a(0.1,1.3); // declaring complex variable a
  Complex b(3.0), c(5.0,-2.3); // dec c vars b and c
  Complex d = b; // dec new c var d
  cout << " d=" << d << ", a=" << a << ", b=" << b << endl;
  d = a*c + b/a; // add, multi and div two c numbers
  cout << "Re(d)=" << d.Re() << ", Im(d)=" << d.Im() << endl; // print real
  // and imaginary parts of d
  return 0;
}
