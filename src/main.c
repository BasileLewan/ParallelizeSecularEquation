#include <omp.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "util.h"

const extern size_t N;

#ifndef N
#define N 10
#endif


// Debug
double RHO = .5;
double DELTA[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
double KSI[10] = {.1, .2, .3, .4, .5, .5, .4, .3, .2, .1};

struct Tuple {
	double a, b, c;
};

double f(double lambda, double rho, double* delta, double* ksi) {
	/* Secular equation, exact value */
	double res = rho;
	for (unsigned short j = 0; j < N; j++)
		res -= sqr(ksi[j]) / (lambda - delta[j]);
	return res;
}

double d_f(double lambda, double* delta, double* ksi) {
	/* Secular equation, first derivative */
	double res = 0;
	for (unsigned short j = 0; j < N; j++) {
		res += (sqr(ksi[j]) / sqr(lambda - delta[j]));
	}
	return res;
}

double d2_f(double lambda, double* delta, double* ksi) {
	/* Secular equation, second derivative */
	double res = 0;
	for (unsigned short j = 0; j < N; j++)
		res -= 2 * sqr(ksi[j]) / pow(lambda - delta[j], 3);
	return res;
}

double d_psi(int k, double lambda, double* delta, double* ksi) {
	double res = 0;
	for (; k > -1; k--)
		res += sqr(ksi[k]) / sqr(lambda - delta[k]);
	return res;
}

struct Tuple middle_way(int k, double lambda, double rho, double* delta, double* ksi, bool fk) {
	/* Parameters to approximate */
	struct Tuple res;
	double fy = f(lambda, rho, delta, ksi);
	double dfy = d_f(lambda, delta, ksi);

	/* interpolating f(x) or fk(x) */
	if (fk) {
		fy += sqr(ksi[k]) / (lambda - delta[k]);
		dfy -= sqr(ksi[k]) / sqr(lambda - delta[k]);
	}

	/* Defining Δ */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = fy - delta_k1 * dfy - d_psi(k, lambda, delta, ksi) * (delta_k - delta_k1);

	return res;
}

struct Tuple fixed_weight(int k, double lambda, double rho, double* delta, double* ksi, bool fk) {
	/* Parameters to approximate */
	struct Tuple res;
	double fy = f(lambda, rho, delta, ksi);
	double dfy = d_f(lambda, delta, ksi);

	if (fk){
		fy += sqr(ksi[k]) / (lambda - delta[k]);
		dfy -= sqr(ksi[k]) / sqr(lambda - delta[k]);
	}

	/* Defining uppercase deltas */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = 0;

	/* Testing if lambda[k] closer to delta[k] or delta[k + 1] */
	if (fabs(lambda - delta[k]) < fabs(lambda - delta[k + 1])) {
		res.c = fy - delta_k1 * dfy - sqr(ksi[k]) * (delta[k] - delta[k + 1]) / sqr(delta_k);
	}
	else{
		res.c = fy - delta_k * dfy - sqr(ksi[k + 1]) * (delta[k + 1] - delta[k]) / sqr(delta_k1);
	}

	return res;
}

struct Tuple gragg(int k, double lambda, double rho, double* delta, double* ksi) {
	/* Defining Δ */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	/* Computing parameters of the Gragg scheme */
	struct Tuple res;
	double fy = f(lambda, rho, delta, ksi);
	double dfy = d_f(lambda, delta, ksi);

	/*** S & s servent à rien ?
	double S = 0, s = 0, c, tmp;
	for (unsigned short i = 0; i < N; i++) {
		if (i == k || i == k + 1)
			continue;
		tmp = sqr(ksi[i]) / pow(delta[i] - lambda, 3);
		s += (delta[i] - delta[k + 1]) * tmp;
		S += (delta[i] - delta[k]) * tmp;
	}
	s *= pow(delta[k] - lambda, 3) / (delta[k] - delta[k + 1]);
	S *= pow(delta[k + 1] - lambda, 3) / (delta[k + 1] - delta[k]);
	s += sqr(ksi[k]); S += sqr(ksi[k + 1]);
	***/

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = f(lambda, rho, delta, ksi) - (delta_k + delta_k1) * d_f(lambda, delta, ksi) + delta_k * delta_k1 * d2_f(lambda, delta, ksi) / 2;

	return res;
}


double initial_guess(int k, double rho, double* delta, double* ksi) {
	/* Defining the first value for the iterations scheme */
	double a, b;
	double fy = f((delta[k] - delta[k + 1]) / 2, rho, delta, ksi);
	double g = fy + 2 * sqr(ksi[k]) / (delta[k + 1] - delta[k]) + 2 * sqr(ksi[k + 1]) / (delta[k] - delta[k + 1]);



	if (k < N - 1) {
		if (fy < 0) {
			a = -g * (delta[k + 1] - delta[k]) + sqr(ksi[k]) + sqr(ksi[k + 1]);
			b = (delta[k] - delta[k + 1]) * sqr(ksi[k + 1]);
			k++;
		}
		else {
			a = g * (delta[k + 1] - delta[k]) + sqr(ksi[k]) + sqr(ksi[k + 1]);
			b = (delta[k + 1] - delta[k]) * sqr(ksi[k]);
		}
		if (a > 0)
			return 2 * b / (a + sqrt(sqr(a) - 4 * b * g)) + delta[k];
		else
			return (a - sqrt(sqr(a) - 4 * b * g)) / (2 * g) + delta[k];
	}
	else {
		double h = fy - g;
		if (fy <= 0 && g <= h) {
			double y = delta[k];
			for (unsigned i = 0; i < N; i++)
				y += sqr(ksi[i]);
			return y;
		}
		else {
			a = g * (delta[k - 1] - delta[k]) + sqr(ksi[k]) + sqr(ksi[k - 1]);
			b = (delta[k - 1] - delta[k]) * sqr(ksi[k]);
			if (a < 0)
				return 2 * b / (a + sqrt(sqr(a) - 4 * b * g)) + delta[k];
			else
				return (a - sqrt(sqr(a) - 4 * b * g)) / (2 * g) + delta[k];
		}
	}
}

double eig_gragg(int k, double rho, double* delta, double* ksi) {
	/* Compute the kth eigenvalue using Gragg's scheme */
	double a, b, c, y, eta = 0;
	struct Tuple param;
	y = initial_guess(k, rho, delta, ksi);

	do {
		y += eta;
		param = gragg(k, y , rho, delta, ksi);
		a = param.a; b = param.b; c = param.c;
		// printf("a = %g b = %g c=%g\n", a, b, c);
		if (a > 0)
			eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c))); // root might be negative for some reason TODO:wtf
		else
			eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
		// printf("y = %g eta = %g\n", y, fabs(eta));
		printf("y = %g eta = %g\n", y, eta);
		printf("1st value : %f, 2nd value : %f, fabs(f(y, rho, delta, ksi)) : %f, fabs(d_f) : %f\n", fabs(f(y, rho, delta, ksi)),  pow(10, -15) + pow(10, -15) * fabs(y) * fabs(d_f(2 * delta[k] - y, delta, ksi)), fabs(d_f(2 * delta[k] - y + eta, delta, ksi)));
	} while (fabs(f(y, rho, delta, ksi)) > pow(10, -15) + pow(10, -15) * fabs(y) * fabs(d_f(2 * delta[k] - y, delta, ksi)));

	return y;
}

double* solve_gragg(double rho, double* delta, double* ksi) {
	/* Compute all the eigenvalues using Gragg's scheme*/
	double* res = malloc(N * sizeof(double));
	if (res == NULL)
		perror("Memory error");
	int k;

#pragma omp parallel for
	for (k = 0; k < N; k++) {
		*(res+k) = eig_gragg(k, rho, delta, ksi);	//check if memory is rightfully allocated
	}
	return res;
}

double eig_hybrid(int k, double rho, double* delta, double* ksi) {
	/* Compute the kth eigenvalue using the Hybrid scheme*/

	/* Initialization */
	struct Tuple param;
	double a, b, c, e, eta, f_nxt;
	double y = initial_guess(k, rho, delta, ksi);
	double f_prv = f(y, rho, delta, ksi);
	double fk = f_prv + sqr(ksi[k]) / (y - delta[k]);  // f(y) without the kth term in the summation

	/* /!\ using 2 or 3 poles /!\ */
	bool use_fk = fk > 0;

	/* getting a new approximation to 1 */
	param = fixed_weight(k, y, rho, delta, ksi, !use_fk);	//modified cause fk > 0 => two poles using f(x)
	a = param.a; b = param.b; c = param.c;
	if (a > 0)
		eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c)));
	else
		eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
	y += eta;

	/* Calculating next value of f(x) */
	f_nxt = f(y, rho, delta, ksi);
	f_prv = f_nxt;
	fk = f_prv + sqr(ksi[k]) / (y - delta[k]);

	/* Checking if we have to switch to the middle way on the 2nd iteration */
	if (f_nxt < 0 && fabs(f_nxt) > 0.1 * fabs(f_prv))
		use_fk = false;

	/* Converging towards zero */
	do {

		//Using fixed_weight or middle_way
		param = use_fk ? fixed_weight(k, y, rho, delta, ksi, false) : middle_way(k, y, rho, delta, ksi, false);
		a = param.a; b = param.b; c = param.c;

		/* Computing η */
		if (a > 0)
			eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c)));
		else
			eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
		y += eta;
		f_nxt = f(y, rho, delta, ksi);

		/* Switch between the schemes*/
		if (f_nxt * f_prv > 0 && fabs(f_nxt) > 0.1 * fabs(f_prv))
			use_fk = !use_fk;


		/* Stopping criterion */
		e = 0;
		for (int j = 0; j <= k; j++) {
			e += (k - j + 6) * fabs(sqr(ksi[j]) / (delta[j] - y));
		}
		for (int j = k; j <= N - 1; j++) {
			e += (j - k + 5) * fabs(sqr(ksi[j]) / (delta[j] - y));
		}
		e += f_nxt;

		/* debug values */
		//printf("y = %g eta = %g\n", y, eta);
		//printf("1st value : %f, 2nd value : %f, e : %f, fabs(y) : %f, fabs(d_f) : %f\n", fabs(f(y, rho, delta, ksi)),  pow(10, -15) * e + pow(10, -15) * fabs(y) * fabs(d_f(2 * delta[k] - y + eta, delta, ksi)), e, fabs(y), fabs(d_f(2 * delta[k] - y + eta, delta, ksi)));

	} while (fabs(f(y, rho, delta, ksi)) > DBL_EPSILON * e + DBL_EPSILON * fabs(y) * fabs(d_f(2 * delta[k] - y + eta, delta, ksi)));

	return y;
}

int main() {
	int K = 5;
	/*
	double* DELTA_100 = malloc(N*sizeof(double));
	double* KSI_100 = malloc(N*sizeof(double));
	double a, sum = 0;

	for (unsigned i = 0; i < N; i++) {
		a = rand();
		DELTA_100[i] = i;
		KSI_100[i] = a;
		sum += sqr(a);
	}
	for (unsigned i = 0; i < N; i++)
		KSI_100[i] /= sum;


	double dn = DELTA_100[N - 1];
	for (unsigned i = 0; i < N; i++)
		dn += KSI_100[i];
	DELTA_100[N] = dn;

	double* tst = solve_gragg(RHO, DELTA_100, KSI_100);

	// double** A = tridiag();

	// print(A);

	// printf("%f\n",tst);
	for (int i = 0; i < N; i++)
		printf("%f\n", tst[i]);
	//free(tst);
	*/

	double y = eig_hybrid(K, RHO, DELTA, KSI);

	return EXIT_SUCCESS;
}
