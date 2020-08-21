// this is a test c++ program

/* this is how to comment in C++*/
using namespace std;
# include <iostream>
# include <cstdlib>
# include <cmath>

int main (int argc, char* argv[])
{
// exception handling, aborts if too few command line args
    if( argc <= 1 ){
      cout << "bad usage:" << argv[0] <<
      " write also a number on the same line, example: prog.exe 0.785" << endl;
      exit(1); // the program exits
    }

// convert the text argv[1] to double using atof:
  double r = atof(argv[1]);
  double s = sin(r);
  cout << "hello, world! sin(" << r <<")=" << s << endl;
// success
  return 0;
}
