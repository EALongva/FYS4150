// header file for the complex class tutorial

#ifndef Complex_H
#define Complex_H
//various include statements and definitions
#include <iostream>
// #include <new>
//#include ....

using namespace std;

class Complex
{
 // definition of variables and their character
private:
  double re, im; // real and imaginary part
public:
  Complex ();                         // Complex c
  Complex (double re, double im = 0.0); // Definition of a complex var
  Complex (const Complex& c);   // usage: Complex(a) // equate two complex vars
  Complex& operator= (const Complex& c); //c=a//equate two c vars, same as prev
  ~Complex () {}                // destructor
  double Re () const; //double real_part = a.Re()
  double Im () const; //double imag_part = a.Im()
  double abs () const;//double m = a.abs() //modulus
  friend Complex operator+ (const Complex& a, const Complex& b);
  friend Complex operator- (const Complex& a, const Complex& b);
  friend Complex operator* (const Complex& a, const Complex& b);
  friend Complex operator/ (const Complex& a, const Complex& b);
};
// declarations of various functions used by the class
#endif
