// function declarations project 2
// https://stackoverflow.com/questions/25274312/is-it-a-good-practice-to-define-c-functions-inside-header-files

#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>

arma::mat maketridiag(double a, double b, double c, int N); // dont need variable names
//a upper off-diagonal, b diagonal and c lower off-diagonal elements, N size

arma::mat BBmatrix(int N);
/*initialize a matrix for the buckling beam with static rho0 and rhoN values,
  rho0 = 0.0 and rhoN = 1.0, N gives the dimensionality of the matrix
*/

arma::mat HOmatrix(double rho0, double rhoN, int N);
/* take rho0 and rhoN, returns an N times N matrix with e = -1/h^2 (h
is the step length) and d = 2/h^2 + V_i on the diagonal. V_i is
computed from a linspace(rho0, rhoN, N) where every element is
squared*/

void maximum_indices(arma::mat A, int N, int& k, int& l);
// finding the maximum value and its indices from matrix A of size N x N

void rotation(arma::mat& A, arma::mat& R, int N, int k, int l);
/*performing a jacobi method rotation on matrix A, overloading the start matrix
  A and the eigenvector matrix R (thus passing them on). The input arguments k
  and l are the indices obtain from the maximum_indices function.*/

arma::mat jacobimethod(arma::mat A, arma::mat& R, int N, int eps, int& iterations);
/*Performs the maximum_indices function and the jacobi rotation in a
  while loop until max = A(k,l) < 10^eps. Returns the diagonalized
  matrix where all off diag elements are below the threshold 10^eps.
  also returns the iterations (how many times the rotation function
  has been called). */

void ToFile(arma::mat A, std::vector<std::string> v, std::string filename);
/*Table of n rows and N columns, array with strings giving
  name of the variables in the table, need to put all values in the
  matrix A so that shape(A) is N x n. e.g v = (dimensionality, iterations,
  time) A is then N x 3 */

void SimTransCount(int eps, int N_min, int N_step, int N_points, std::string filename);
/*FOR THE BUCKLING BEAM MATRIX, obtained with 'maketridiag'
  Runs the jacobimethod for linearly spaced values of N (dimensionality of matrix)
  and counts the iterations of similarity transformations needed to get the off -
  diagonal elements below the threshold. Also calculate the CPU time for each
  jacobimethod run. Saves the results to a file, name given in filename (using ToFile).
*/


#endif
