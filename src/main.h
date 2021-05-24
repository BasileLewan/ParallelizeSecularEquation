#pragma once
#include <stdbool.h>

/* Secular function and additional functions */
double f(double lambda, double rho, double* delta, double* ksi);
double d_f(double lambda, double* delta, double* ksi);
double fk(int k, double lambda, double rho, double* delta, double* ksi);
double d_fk(int k, double lambda, double* delta, double* ksi);
double d2_f(double lambda, double* delta, double* ksi);
double d_psi(int k, double lambda, double* delta, double* ksi);
double g(int k, double lambda, double rho, double* delta, double* ksi);

/* Gragg algorithm */
struct Tuple gragg(int k, double lambda, double rho, double* delta, double* ksi);
double eig_gragg(int k, double rho, double* delta, double* ksi);
double* solve_gragg(double rho, double* delta, double* ksi);

/* Middle-way algorithm */
struct Tuple middle_way(int k, double lambda, double rho, double* delta, double* ksi, bool pol2);
struct Tuple fixed_weight(int k, double lambda, double rho, double* delta, double* ksi, bool pol2);
double initial_guess(int k, double rho, double* delta, double* ksi);
double eig_hybrid(int k, double rho, double* delta, double* ksi, int* compteur, double precision);
double* solve_hybrid(double rho, double* delta, double* ksi, int* compteur, double precision);


