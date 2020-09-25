// function declarations project 2
// https://stackoverflow.com/questions/25274312/is-it-a-good-practice-to-define-c-functions-inside-header-files

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>

arma::mat maketridiag(double a, double b, double c, int N); // dont need variable names
//a upper off-diagonal, b diagonal and c lower off-diagonal elements, N size

void maximum_indices(arma::mat A, int N, int& k, int& l);
// finding the maximum value and its indices from matrix A of size N x N

arma::mat rotation(arma::mat A, int N, int k, int l);
// performing a jacobi method rotation on matrix A, returning the changed
// matrix A. The input arguments k and l are the indices obtain from the
// maximum_indices function.

arma::mat jacobimethod(arma::mat A, int N, int eps, int& iterations);
// Performs the maximum_indices function and the jacobi rotation in a
// while loop until max = A(k,l) < 10^eps. Returns the diagonalized
// matrix where all off diag elements are below the threshold 10^eps.
// also returns the iterations (how many times the rotation function
// has been called).

void ToFile(arma::mat A, arma::vec v, int N, char filename);
// maybe: table of n rows and N columns, array with strings giving
// name of the variables in the table, need to put all values in the
// matrix A so that shape(A) is N x n. e.g v = (dimensionality, iterations,
// time) A is then N x 3


#endif
