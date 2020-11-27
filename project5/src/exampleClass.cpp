#include "../include/exampleClass.hpp"

using namespace arma;
using namespace std;

exampleClass::exampleClass(int N)
{
    number = N;
    std::cout << "init number = " << number << std::endl;
}

void exampleClass::function()
{
  std::cout << "test function for exampleClass, number = " << number << std::endl;
}

//
