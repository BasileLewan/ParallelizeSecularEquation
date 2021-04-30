#include <omp.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "util.h"

/* N defines the matices length */
const extern size_t N;

#ifndef N
#define N 50
#endif


/* Initialization */
double RHO = .5;
double DELTA[N + 1] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0};
double KSI[N] = {.1, .2, .3, .4, .5, .5, .4, .3, .2, .1};

/* A structure for the 3 parameters we use to interpolate f(x) */
struct Tuple {
	double a, b, c;
};

/* Secular equation, exact value */
double f(double lambda, double rho, double* delta, double* ksi) {
	double res = rho;
	for (unsigned short j = 0; j < N; j++)
		res -= sqr(ksi[j]) / (lambda - delta[j]);
	return res;
}

/* Secular equation, first derivative */
double d_f(double lambda, double* delta, double* ksi) {
	double res = 0;
	for (unsigned short j = 0; j < N; j++) {
		res += (sqr(ksi[j]) / sqr(lambda - delta[j]));
	}
	return res;
}

/* Secular equation without the kth term in the summation */
double fk(int k, double lambda, double rho, double* delta, double* ksi) {
	double res = rho;
	for (unsigned short j = 0; j < N; j++) {
		if (j == k)
			continue;
		res -= sqr(ksi[j]) / (lambda - delta[j]);
			}
	return res;
}

/* Secular equation without the kth term in the summation, first derivative */
double d_fk(int k,  double lambda, double* delta, double* ksi) {
	double res = 0;
	for (unsigned short j = 0; j < N; j++) {
		if (j == k)
			continue;
		res += (sqr(ksi[j]) / sqr(lambda - delta[j]));
	}
	return res;
}

/* Secular equation, second derivative */
double d2_f(double lambda, double* delta, double* ksi) {
	double res = 0;
	for (unsigned short j = 0; j < N; j++)
		res -= 2 * sqr(ksi[j]) / pow(lambda - delta[j], 3);
	return res;
}

/* Function Psi, first derivative
	 Used to calculate parameter c in the middle way method */
double d_psi(int k, double lambda, double* delta, double* ksi) {
	double res = 0;
	for (; k > -1; k--)
		res += sqr(ksi[k]) / sqr(lambda - delta[k]);
	return res;
}

/* The middle way method
   Solves the secular equation using both nearby poles */
struct Tuple middle_way(int k, double lambda, double rho, double* delta, double* ksi, bool pol2) {
	/* Initialization */
	struct Tuple res;
	double fy = pol2 ? fk(k,lambda, rho, delta, ksi) :f(lambda, rho, delta, ksi);
	double dfy = pol2 ? d_fk(k, lambda, delta, ksi) :d_f(lambda, delta, ksi);

	/* Defining Δ */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = fy - delta_k1 * dfy - d_psi(k, lambda, delta, ksi) * (delta[k] - delta[k + 1]);

	return res;
}

/* The fixed weight method
   Solves the secular equation using the most appropriated pole */
struct Tuple fixed_weight(int k, double lambda, double rho, double* delta, double* ksi, bool pol2) {
	/* Initialization */
	struct Tuple res;
	double fy = pol2 ? fk(k,lambda, rho, delta, ksi) :f(lambda, rho, delta, ksi);
	double dfy = pol2 ? d_fk(k, lambda, delta, ksi) :d_f(lambda, delta, ksi);

	/* Defining Δ */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = 0;

	/* Testing if lambda[k] closer to delta[k] or delta[k + 1]
		 Then uses most appropriated pole */
	if (fabs(lambda - delta[k]) < fabs(lambda - delta[k + 1])) {
		res.c = fy - delta_k1 * dfy - sqr(ksi[k]) * (delta[k] - delta[k + 1]) / sqr(delta_k);
	}
	else {
		res.c = fy - delta_k * dfy - sqr(ksi[k + 1]) * (delta[k + 1] - delta[k]) / sqr(delta_k1);
	}

	return res;
}

/* Gragg's Scheme
	 Uses f second derivative to calculate parameter c */
struct Tuple gragg(int k, double lambda, double rho, double* delta, double* ksi) {
	/* Initialization */
	struct Tuple res;
	double fy = f(lambda, rho, delta, ksi);
	double dfy = d_f(lambda, delta, ksi);

	/* Defining Δ */
	double delta_k = delta[k] - lambda;
	double delta_k1 = delta[k + 1] - lambda;

	res.a = (delta_k + delta_k1) * fy - delta_k * delta_k1 * dfy;
	res.b = delta_k * delta_k1 * fy;
	res.c = f(lambda, rho, delta, ksi) - (delta_k + delta_k1) * d_f(lambda, delta, ksi) + delta_k * delta_k1 * d2_f(lambda, delta, ksi) / 2;

	return res;
}

/* Defining the first value for the iterations scheme */
double initial_guess(int k, double rho, double* delta, double* ksi) {
	/* Initialization */
	double a, b;

	/* f((dk + dk+1) / 2) */
	double fy = f((delta[k] + delta[k + 1]) / 2, rho, delta, ksi);

	/* g((dk + dk+1) / 2), with f(x) = g(x) + h(x) */
	double g = fy + 2 * sqr(ksi[k]) / (delta[k + 1] - delta[k]) + 2 * sqr(ksi[k + 1]) / (delta[k] - delta[k + 1]);

	if (k < N - 1) {
		/* 2 poles available */
		/* calculating parameters a and b */
		if (fy < 0) {
			a = -g * (delta[k + 1] - delta[k]) + sqr(ksi[k]) + sqr(ksi[k + 1]);
			b = (delta[k + 1] - delta[k]) * sqr(ksi[k + 1]);
			k++;
		}
		else {
			a = g * (delta[k + 1] - delta[k]) + sqr(ksi[k]) + sqr(ksi[k + 1]);
			b = (delta[k + 1] - delta[k]) * sqr(ksi[k]);
		}
		if (a > 0) {
			return (2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * g)))) + delta[k]; // negative sqrt !!
		}
		else {
			return ((a - sqrt(fabs(sqr(a) - 4 * b * g))) / (2 * g)) + delta[k]; // negative sqrt !!
		}
	}
	else {
		/* last interval */

		/* defining dn+1 as dn + ztz/rho as in the paper */
		double dn = 0;
		for (int i = 0; i < N; ++i)
			dn += sqr(ksi[i]);
		dn /= rho;
		dn += delta[k];

		/* f((dn + dn+1) / 2) */
		fy = f((delta[k] + dn) / 2, rho, delta, ksi);

		/* g((dn + dn+1) / 2), with f(x) = g(x) + h(x) */
		g = rho;
		for (int i = 0; i < N - 1; ++i)
			g -= sqr(ksi[i]) / ((delta[k] + dn) / 2 - delta[i]);

		/* f(dn) */
		double fn1 = f(dn, rho, delta, ksi);

		/* g(dn) */
		double gn1 = rho;
		for (int i = 0; i < N - 1; ++i)
			gn1 -= sqr(ksi[i]) / (dn - delta[i]);

		/* h(dn) */
		double h = fn1 - gn1;

		/* g((dn-1 + dn) / 2) */
		double g_1 = fy + 2 * sqr(ksi[k - 1]) / (delta[k] - delta[k - 1]) + 2 * sqr(ksi[k]) / (delta[k - 1] - delta[k]);

		if (fy <= 0 && g_1 <= -h) {
			//printf("case 1\n");
			return dn;
		}
		else {
			//printf("case 2\n");
			a = g * (delta[k - 1] - delta[k]) + sqr(ksi[k]) + sqr(ksi[k - 1]);
			b = (delta[k - 1] - delta[k]) * sqr(ksi[k]);
			if (a < 0)
				return 2 * b / (a + sqrt(sqr(a) - 4 * b * g)) + delta[k];
			else
				return (a - sqrt(sqr(a) - 4 * b * g)) / (2 * g) + delta[k];
		}
	}
}


/* Compute the kth eigenvalue using Gragg's scheme */
double eig_gragg(int k, double rho, double* delta, double* ksi) {
	/* Initialization */
	double a, b, c, y, eta = 0;
	struct Tuple param;
	y = initial_guess(k, rho, delta, ksi);

	do {
		y += eta;
		param = gragg(k, y , rho, delta, ksi);
		a = param.a; b = param.b; c = param.c;
		if (a > 0)
				eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c)));
		else
			eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
	} while (fabs(f(y, rho, delta, ksi)) > pow(10, -15) + pow(10, -15) * fabs(y) * fabs(d_f(2 * delta[k] - y, delta, ksi)));

	return y;
}


/* Compute all the eigenvalues using Gragg's scheme*/
double* solve_gragg(double rho, double* delta, double* ksi) {
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


/* Compute the kth eigenvalue using the Hybrid scheme*/
double eig_hybrid(int k, double rho, double* delta, double* ksi, int* compteur, double precision) {
	/* Initialization */
	struct Tuple param;
	double a, b, c, e, eta, f_nxt;
	double y = initial_guess(k, rho, delta, ksi);
	double ig = y;
	double f_prv = f(y, rho, delta, ksi);
	double f_k = fk(k, y, rho, delta, ksi);  // f(y) without the kth term in the summation
	bool swtch = true;

	/* /!\ using 2 or 3 poles /!\ */
	bool use_fk = f_k < 0;

	/* getting a new approximation to 1 */
	param = fixed_weight(k, y, rho, delta, ksi, use_fk);	//modified cause fk > 0 => two poles using f(x)
	a = param.a; b = param.b; c = param.c;
	if (a > 0)
		eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c)));
	else
		eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
	y += eta;

	/* Calculating next value of f(x) */
	f_nxt = f(y, rho, delta, ksi);
	// f_k = f_prv + sqr(ksi[k]) / (y - delta[k]);

	/* Checking if we have to switch to the middle way on the 2nd iteration */
	if (f_nxt < 0 && fabs(f_nxt) > 0.1 * fabs(f_prv))
		swtch = false;

	/* Testing if lambda[k] closer to delta[k] or delta[k + 1] */
	int kj = fabs(y - delta[k]) < fabs(y - delta[k + 1]) ? k : k + 1;

	/* debug values */
	int i = 1;
	int count_switch = 0;

	/* Converging towards zero */
	do {
		if (count_switch > 5) {
			/* more than 5 scheme switches in a row */
			/* breaking and returning current y */

			/* debug values for infinite iterations */
			/*
			printf("\n------------------------- k = %.2d -------------------------\n\n", k);
			printf("delta[k-1] : %.1f, delta[k] : %.1f, delta[k+1] : %.1f\n", delta[k - 1], delta[k], delta[k + 1]);
			printf("ksi[k-1] : %.3f, ksi[k] : %.3f, ksi[k+1] : %.3f\n", ksi[k - 1], ksi[k], ksi[k + 1]);
			printf("initial guess : %.2f\n", ig);
			printf("%e > %e\n", fabs(f(y, rho, delta, ksi)), precision * e + precision * fabs(y - delta[kj]) * fabs(d_f(2 * delta[kj] - y, delta, ksi)));
			printf("y : %f, f(y) : %e\n", y, f(y, rho, delta, ksi));
			printf("\n----------------------------------------------------------\n\n");
			*/
			*compteur += 1;
			break;
		}

		/* Using fixed_weight or middle_way */
		param = swtch ? fixed_weight(k, y, rho, delta, ksi, false) : middle_way(k, y, rho, delta, ksi, false);
		a = param.a; b = param.b; c = param.c;

		/* Computing η */
		if (a > 0)
			eta = 2 * b / (a + sqrt(fabs(sqr(a) - 4 * b * c)));
		else
			eta = (a - sqrt(fabs(sqr(a) - 4 * b * c))) / (2 * c);
		y += eta;
		f_prv = f_nxt;
		f_nxt = f(y, rho, delta, ksi);

		/* Switch between the schemes*/
		if (f_nxt * f_prv > 0 && fabs(f_nxt) > 0.1 * fabs(f_prv)) {
			swtch = !swtch;
			count_switch ++;
		}else{
			count_switch = 0;
		}

		/* Stopping criterion */
		e = 2 * rho;
		for (int j = 0; j <= k; j++) {
			e += (k - j + 6) * fabs(sqr(ksi[j]) / (delta[j] - y));
		}
		for (int j = k + 1; j <= N - 1; j++) {
			e += (j - k + 5) * fabs(sqr(ksi[j]) / (delta[j] - y));
		}
		e += f_prv;

		++i;

	} while (fabs(f(y, rho, delta, ksi)) > precision * e + precision * fabs(y - delta[kj]) * fabs(d_f(2 * delta[kj] - y, delta, ksi)));

	return y;

}


/* Compute all the eigenvalues using the Hybrid scheme*/
double* solve_hybrid(double rho, double* delta, double* ksi, int* compteur, double precision) {
	double* res = malloc(N * sizeof(double));
	if (res == NULL)
		perror("Memory error");
	int k;

#pragma omp parallel for
	for (k = 0; k < N; k++) {
		*(res+k) = eig_hybrid(k, rho, delta, ksi, compteur, precision);	//check if memory is rightfully allocated
	}
	return res;
}


/* Tests various precisions on matrices of various length
	.pre0    : precision to start with and then lower
	.preInc  : number of precision incrementation before lowering to next exponential
		ex : if one choses pre0 = 1e-15 and preInc = 3
			 	 first 3 iterations would have a precision equal to 1e-15, 3.3e-15 and 6.6e-15
				 fourth iteration would get to a precision equal to 1e-14, etc.
	.nbItPre : number of iterations on precision before the method stops
	.len0    : length of matrices to start with and then increase
	.lenInc  : length to add to matrices at each iteration
	.nbItLen : number of iterations on matrices length before lowering precision
	.nbTests : number of tests on matrices on given length with given precision

	test_precision first increases matrices length nbItLen times for a given precision,
	then lowers precision and reset matrices length, and repeats this process preInc times */
void test_precision(double pre0, double preInc, int nbItPre, int len0, int lenInc, int nbItLen, int nbTests) {
	double pre = pre0;
	int len = len0;
	int quo = (int) preInc;
	int compteur = 0;
	int compteur1 = 0;
	int *cpt = &compteur;
	int nbTestsPerPre = 0;
	for (int i = 0; i < nbItLen; ++i)
		nbTestsPerPre += nbTests * (len0 + i * lenInc);
	double* DELTA_LEN = malloc(len*sizeof(double));
	double* KSI_LEN = malloc(len*sizeof(double));
	double RHO = .5;
	double ksisum = 0;
	srand(time(NULL));

	for (int cpt1 = 0; cpt1 < nbItPre; ++cpt1) {
		/* It on precision */
		compteur1 = 0;
		if (cpt1 % quo == 0)
			pre = pre0 * pow(10, cpt1 / quo);
		else
			pre = pre0 * (10 / preInc) * (cpt1 % quo) * pow(10, cpt1 / quo);

		printf("/---------------------------------------------------------------------/\n");
		printf("/------------------------ precision : %.1e ------------------------/\n", pre);
		printf("/---------------------------------------------------------------------/\n");

		for (int cpt2 = 0; cpt2 < nbItLen; ++cpt2) {
			/* It on matrices length */
			compteur = 0;
			len = len0 + cpt2 * lenInc;

			printf("         --------------- matrices length : %.2d ----------------\n", len);

			for (int cpt3 = 0; cpt3 < nbTests; ++cpt3) {
				/* Testing precision pre on nbTests generated matrices length len */
				/* Initialization */
				ksisum = 0;

				/* Generating DELTA and KSI matrices */
				for (unsigned i = 0; i < len; i++) {
					if(i == 0)
						KSI_LEN[i] = rand() % 10;
					else
						KSI_LEN[i] = DELTA_LEN[i - 1] + (rand() % 5) + 1;	//delta <= 100 for len <= 50
					KSI_LEN[i] = ((rand() % 95) + 5) * pow(10, -2);				//ksi is in ]0; 1[
					ksisum += sqr(KSI_LEN[i]);
				}

				/* Defining d[n] and ksi[n] for last iteration */
				KSI_LEN[len] = ((rand() % 95) + 5) * pow(10, -2);
				DELTA_LEN[len] = DELTA_LEN[len - 1] + ksisum / RHO;	//defining dn+1 as in the paper : dn + ztz / rho

				solve_hybrid(RHO, DELTA_LEN, KSI_LEN, cpt, pre);
			}
			compteur1 += compteur;
			printf("                             breaks : %.2d/%.2d                   \n", compteur, len * nbTests);
		}
		printf("                        total breaks : %.2d/%.2d                   \n", compteur1, nbTestsPerPre);
	}
}


int main() {
	/* test Lewan sur matrices */

	test_precision(pow(10, -15), 2, 8, 40, 20, 6, 3);

	/*
	//Initialization
	double* DELTA_100 = malloc(N*sizeof(double));
	double* KSI_100 = malloc(N*sizeof(double));
	double ksisum = 0;
	srand(time(NULL));

	//Generating DELTA and KSI matrices
	for (unsigned i = 0; i < N; i++) {
		if(i == 0)
			DELTA_100[i] = rand() % 10;
		else
			DELTA_100[i] = DELTA_100[i - 1] + (rand() % 5) + 1;	//delta <= 100 for N <= 50
		KSI_100[i] = ((rand() % 95) + 5) * pow(10, -2);				//ksi is in ]0; 1[
		ksisum += sqr(KSI_100[i]);
	}

	//Defining d[n] and ksi[n] for last iteration
	KSI_100[N] = ((rand() % 95) + 5) * pow(10, -2);
	DELTA_100[N] = DELTA_100[N - 1] + ksisum / RHO;	//defining dn+1 as in the paper : dn + ztz / rho

	//Counting number of algorithm forced stops
	int compteur = 0;
	int *cpt = &compteur;

	//Precision used in stopping criteria
	double precision = 5 * pow(10, -15);

	double *hyb = solve_hybrid(RHO, DELTA_100, KSI_100, cpt, precision);
	printf("          nb de break : %d\n\n", compteur);

	printf("----- précision utilisée : %.1e -----\n", precision);
	*/

	return EXIT_SUCCESS;
}
