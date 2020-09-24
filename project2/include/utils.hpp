// function declarations project 2
// https://stackoverflow.com/questions/25274312/is-it-a-good-practice-to-define-c-functions-inside-header-files

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <armadillo>

arma::mat maketridiag(double a, double b, double c, int N); // dont need variable names
// b diagonal and a, c off diagonal elements, N size

void maximum_indices(arma::mat A, int N, int& k, int& l);
// finding the maximum value and its indices from matrix A of size N x N


#endif
