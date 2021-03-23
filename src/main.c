#include <math.h>
#include <float.h>
#include <stdio.h>
// #include <unistd.h>
#define D 10

const double RHO = .5;
const double DELTA[D] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
const double KSI[D] = { .1, .2, .3, .4, .5, .5, .4, .3, .2, .1 };

double sqr(double x) {
	return x * x;
}
struct Tuple {
	double a, b, c;
};

double f(double lambda, double rho, double* delta, double* ksi) {
	/* Secular equation, exact value */
	double res = rho;
	for (unsigned short j = 0; j < D; j++)
		res -= sqr(ksi[j]) / (lambda - delta[j]);
	return res;
}

double d_f(double lambda, double rho, double* delta, double* ksi) {
	/* Secular equation, first derivative */
	double res = 0;
	for (unsigned short j = 0; j < D; j++)
		res += sqr(ksi[j]) / sqr(lambda - delta[j]);
	return res;
}

double d2_f(double lambda, double rho, double* delta, double* ksi) {
	/* Secular equation, second derivative */
	double res = 0;
	for (unsigned short j = 0; j < D; j++)
		res -= 2 * sqr(ksi[j]) / pow(lambda - delta[j], 3);
	return res;
}

double d_psi(int k, double lambda, double rho, double* delta, double* ksi) {
	double res = 0;
	for (; k > -1; k--)
		res += sqr(ksi[k]) / sqr(lambda - delta[k]);
	return res;
}

struct Tuple middle_way(int k, double lambda, double rho, double* delta, double* ksi) {
	/* Parameters to approximate */
	struct Tuple res;
	double fy = f(lambda, rho, delta, ksi);
	double dfy = d_f(lambda, rho, delta, ksi);

	/* Defining uppercase deltas */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = fy - delta_k1 * dfy - d_psi(k, lambda, rho, delta, ksi) * (delta_k - delta_k1);

	return res;
}

struct Tuple fixed_weight(int k, double lambda, double rho, double* delta, double* ksi) {
	/* Parameters to approximate */
	struct Tuple res;
	double fy = f(lambda, rho, delta, ksi);
	double dfy = d_f(lambda, rho, delta, ksi);

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
	/* Defining uppercase deltas */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	/* Computing parameters of the Gragg scheme */
	double S = 0, s = 0, c, tmp;
	for (unsigned short i = 0; i < D; i++) {
		if (i == k || i == k + 1)
			continue;
		tmp = sqr(ksi[i]) / pow(delta[i] - lambda, 3);
		s += (delta[i] - delta[k + 1]) * tmp;
		S += (delta[i] - delta[k]) * tmp;
	}
	s *= pow(delta[k] - lambda, 3) / (delta[k] - delta[k + 1]);
	S *= pow(delta[k + 1] - lambda, 3) / (delta[k + 1] - delta[k]);
	s += sqr(ksi[k]); S += sqr(ksi[k + 1]);
	c = f(lambda, rho, delta, ksi) - (delta_k + delta_k1) * d_f(lambda, rho, delta, ksi) + delta_k * delta_k1 * d2_f(lambda, rho, delta, ksi) / 2;
	struct Tuple res;
	res.a = s; res.b = S; res.c = c;
	return res;
}

double initial_guess(int k, double rho, double* delta, double* ksi) {
	/* Defining the first value for the iterations scheme */
	// TODO : case >λₙ
	double a, b;
	double fy = f((delta[k] - delta[k + 1]) / 2, rho, delta, ksi);
	double g = fy +  2 * sqr(ksi[k]) / (delta[k + 1] - delta[k]) + 2 * sqr(ksi[k + 1]) / (delta[k] - delta[k + 1]);
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
		return 2 * b / (a + sqrt(sqr(a) - 4 * b * g));// + delta[k]
	else
		return (a - sqrt(sqr(a) - 4 * b * g)) / (2 * g);// + delta[k]
}

int main() {
	struct Tuple param, grg;
	double a, b, c, eta = 0;
	int K = 4;
	double y = initial_guess(K, RHO, DELTA, KSI);

	do {
		y += eta;
		param = middle_way(K, y, RHO, DELTA, KSI);
		grg = gragg(K, y, RHO, DELTA, KSI);
		a = param.a; b = param.b; c = grg.c;
		printf("a = %g b = %g c=%g\n",a, b, c);
		if (a > 0)
			eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c))); // root might be negative for some reason TODO:wtf
		else
			eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
		printf("y = %g eta = %g\n", y, fabs(eta));
		//printf("stop = %g\n", 16. * FLT_MIN * fmin(fabs(DELTA[K] - y), fabs(DELTA[K + 1] - y)));
		//printf("y = %g eta = %g f(delta + to) = %g, stop = %g\n", y, fabs(eta), fabs(f(y, RHO, DELTA, KSI)), pow(10, -15) + pow(10, -15) * fabs(y) * fabs(d_f(2 * DELTA[K] - y, RHO, DELTA, KSI)));
	} while (fabs(f(y, RHO, DELTA, KSI)) > pow(10, -15) + pow(10, -15) * fabs(y) * fabs(d_f(2 * DELTA[K] - y, RHO, DELTA, KSI)));
		//fabs(eta) > 16. * FLT_MIN * fmin(fabs(DELTA[K] - y), fabs(DELTA[K + 1] - y)) // naive stopping criterion

	return 1;
}
