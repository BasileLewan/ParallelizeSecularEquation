#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include "main.h"
#include "util.h"
extern size_t N;

#define GRAGG 0
#define HYBRID 1

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
	unsigned len = len0;
	int quo = (int)preInc;
	int compteur = 0;
	int compteur1 = 0;
	int* cpt = &compteur;
	int nbTestsPerPre = 0;
	for (int i = 0; i < nbItLen; ++i)
		nbTestsPerPre += nbTests * (len0 + i * lenInc);
	double* DELTA_LEN = malloc(len * sizeof(double));
	double* KSI_LEN = malloc(len * sizeof(double));
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

				/* Generating DELTA and KSI vectors */
				for (unsigned i = 0; i < len; i++) {
					if (i == 0)
						KSI_LEN[i] = rand() % 10;
					else
						KSI_LEN[i] = DELTA_LEN[i - 1] + (rand() % 5) + 1;	//delta <= 100 for len <= 50
					KSI_LEN[i] = ((rand() % 95) + 5) * pow(10, -2);				//ksi is in ]0; 1[
					ksisum += sqr(KSI_LEN[i]);
				}

				/* Defining d[n] and ksi[n] for last iteration */
				KSI_LEN[len] = ((rand() % 95) + 5) * pow(10, -2);
				DELTA_LEN[len] = DELTA_LEN[len - 1] + ksisum / RHO;	//defining dn+1 as in the paper : dn + ztz / rho

				double* res = solve_hybrid(RHO, DELTA_LEN, KSI_LEN, cpt, pre);
				free(res);
			}
			compteur1 += compteur;
			printf("                             breaks : %.2d/%.2d                   \n", compteur, len * nbTests);
		}
		printf("                        total breaks : %.2d/%.2d                   \n", compteur1, nbTestsPerPre);
	}
}

double speed_test(double** A, unsigned algo) {
	/* running a speedtest on a matrix*/
	double ρ = get_rho(A);
	double* ξ = get_ksi(A);
	double* δ = get_delta(A, ξ, ρ);
	int ct;
	double* eig;
	double start = omp_get_wtime();
	if (algo)
		eig = solve_hybrid(ρ, δ, ξ, &ct, DBL_EPSILON);
	else
		eig = solve_gragg(ρ, δ, ξ);
	return omp_get_wtime() - start;
}


int main() {
	// test_precision(pow(10, -15), 2, 8, 40, 20, 6, 3);

	fprintf(stderr, "Number of available threads : %d\n\n", omp_get_max_threads());

	for (N = 1000; N < 50001; N += 1000) {
		double** A = gen_sym(N);
		fprintf(stderr, "Size : %d\n\tGragg  : %.10f\n", N, speed_test(A, GRAGG));
		fprintf(stderr, "\tHybrid : %.10f\n", speed_test(A, HYBRID);
		free(A[0]);
	}

	return EXIT_SUCCESS;
}