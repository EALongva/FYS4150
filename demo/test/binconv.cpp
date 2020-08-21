// program that converts int number into binary

using namespace std;
# include <iostream>
# include <cstdlib>
# include <cmath>

int main (int argc, char* argv[])
{
  int i;
  int terms[32]; // list that stores the binary numbers
  int number = atoi(argv[1]);
// initialize or compute the terms a0, a1 etc using a for loop
  for (i=0; i < 32 ; i++){
    terms[i] = 0;
  }
  for (i=0; i < 32 ; i++){
    terms[i] = number%2;
    number /= 2;
  }
// write out results
  cout << " Number of bytes used= " << sizeof(number) << endl;
  for (i=0; i < 32 ; i++){
    cout << " Term nr. " << i << " Value= " << terms[i];
    cout << terms[::] << endl;
  }
  return 0;
}
