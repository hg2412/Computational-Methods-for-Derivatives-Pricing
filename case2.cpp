/*
 * case2.cpp

 *
 *  Created on: Oct 22, 2015
 *      Author: Haoxiang
 */

#define Wat(i,j) (gsl_matrix_get(wMatrix,i,j))

#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <utility>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_expint.h>

using namespace std;

//print matrix
void printmatrix(gsl_matrix * mat, int N, int M) {
	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= M; j++)
			cout << gsl_matrix_get(mat, i, j) << " ";
		cout << endl;
	}

}
//output matrix as .csv file
void writematrix(gsl_matrix * mat, int N, int M) {
	ofstream fout("mat.csv");
	for (int i = 0; i <= N; i++) {
		for (int j = 0; j <= M; j++)
			fout << gsl_matrix_get(mat, i, j) << ",";
		fout << "\n";
	}
	fout.close();
}

//print vector
void printvector(vector<double> v, int N) {
	for (int i = 0; i <= N; i++) {
		cout << v[i] << ",";
	}
	cout << endl;
}

class BarrierOption {
public:
	double S; //spot price
	double K; //strike price
	double B; //upper barrier
	double r; //risk-free rate
	double q; //dividend rate
	double T; // maturity
	double sigma; //volatility
	double vega;
	double theta;
	double Y;
	BarrierOption() {
		S = 1900;
		K = 2000;
		B = 2200;
		r = 0.0025;
		q = 0.015;
		T = 0.5;
		sigma = 0.25;
		vega = 0.31;
		theta = -0.25;
		Y = 0.4;
	}
};

class PIDEsolver: public BarrierOption {
public:
	int M;
	int N;
	double xmax;
	double xmin;
	double tmax;
	double dx;
	double dt;

	gsl_vector * xMesh;
	gsl_vector * tMesh;
	gsl_matrix * wMatrix;

	vector<double> g1klamdap;
	vector<double> g1klamdan;
	vector<double> g2klamdap;
	vector<double> g2klamdan;
	vector<double> g2klamdan_plus1;
	vector<double> g2klamdap_minus1;
	vector<double> g2klamdap_plus1;

	double lamdap() {
		return sqrt(pow(theta, 2) / pow(sigma, 4) + 2 / (vega * pow(sigma, 2)))
				- theta / pow(sigma, 2);
	}

	double lamdan() {
		return sqrt(pow(theta, 2) / pow(sigma, 4) + 2 / (vega * pow(sigma, 2)))
				+ theta / pow(sigma, 2);
	}

	double g1(double alpha, double x) {
		if (alpha == 0)
			return exp(-x);
		if (x == 0)
			return gsl_sf_gamma(1 - alpha);
		return (gsl_sf_gamma_inc(1 - alpha, x));
	}

	double g2(double alpha, double x) {
		if (alpha == 0)
			return gsl_sf_expint_E1(x);
		return exp(-x) * pow(x, -alpha) / alpha
				- gsl_sf_gamma_inc(1 - alpha, x) / alpha;

	}

	double sigma2eps(double eps) {
		double t1 = (1 / vega) * pow(lamdap(), Y - 2)
				* (-pow(lamdap() * eps, 1 - Y) * exp(-lamdap() * eps)
						+ (1 - Y) * (g1(Y, 0) - g1(Y, lamdap() * eps)));
		double t2 = (1 / vega) * pow(lamdan(), Y - 2)
				* (-pow(lamdan() * eps, 1 - Y) * exp(-lamdan() * eps)
						+ (1 - Y) * (g1(Y, 0) - g1(Y, lamdan() * eps)));
		return t1 + t2;

	}

	double weps(double eps) {
		return pow(lamdap(), Y) / vega * g2(Y, lamdap() * eps)
				- pow(lamdap() - 1, Y) / vega * g2(Y, (lamdap() - 1) * eps)
				+ pow(lamdan(), Y) / vega * g2(Y, lamdan() * eps)
				- pow(lamdan() + 1, Y) / vega * g2(Y, (lamdan() + 1) * eps);
	}

	double l() {
		double Bl = sigma2eps(dx) * dt / (2 * dx * dx)
				- (r - q + weps(dx) - 0.5 * sigma2eps(dx)) * dt / (2 * dx);
		return -Bl;
	}

	double u() {
		double Bu = sigma2eps(dx) * dt / (2 * dx * dx)
				+ (r - q + weps(dx) - 0.5 * sigma2eps(dx)) * dt / (2 * dx);
		return -Bu;
	}

	double d(int N, int i) {
		double Bl = sigma2eps(dx) * dt / (2 * dx * dx)
				- (r - q + weps(dx) - 0.5 * sigma2eps(dx)) * dt / (2 * dx);
		double Bu = sigma2eps(dx) * dt / (2 * dx * dx)
				+ (r - q + weps(dx) - 0.5 * sigma2eps(dx)) * dt / (2 * dx);
		return 1 + r * dt + Bl + Bu
				+ dt / vega
						* (pow(lamdan(), Y) * g2(Y, i * dx * lamdan())
								+ pow(lamdap(), Y)
										* g2(Y, (N - i) * dx * lamdap()));
	}

	double R(int i, int j) {
		double sum = 0;
		for (int k = 1; k <= i - 1; k++) {
			sum = sum
					+ pow(lamdan(), Y)
							* (Wat(i - k, j) - Wat(i, j)
									- k * (Wat(i-k-1,j) - Wat(i - k, j)))
							* (g2klamdan[k] - g2klamdan[k + 1]);
		}
		for (int k = 1; k <= i - 1; k++) {
			sum = sum
					+ (Wat(i-k-1,j) - Wat(i - k, j))
							/ (pow(lamdan(), 1 - Y) * dx)
							* (g1klamdan[k] - g1klamdan[k + 1]);
		}
		for (int k = 1; k <= N - i - 1; k++) {
			sum = sum
					+ pow(lamdap(), Y)
							* (Wat(i + k, j) - Wat(i, j)
									- k * (Wat(i+k+1,j) - Wat(i + k, j)))
							* (g2klamdap[k] - g2klamdap[k + 1]);
		}

		for (int k = 1; k <= N - i - 1; k++) {
			sum = sum
					+ (Wat(i+k+1,j) - Wat(i + k, j))
							/ (pow(lamdap(), 1 - Y) * dx)
							* (g1klamdap[k] - g1klamdap[k + 1]);
		}
		return sum;
	}

	PIDEsolver(int m, int n) {
		BarrierOption();
		M = m;
		N = n;
		xmin = log(1000);
		xmax = log(B);
		dx = (xmax - xmin) / N;
		tmax = T;
		dt = tmax / M;
		xMesh = gsl_vector_alloc(N + 1);
		tMesh = gsl_vector_alloc(M + 1);
		//Initialize SMeth
		for (int i = 0; i <= N; i++)
			gsl_vector_set(xMesh, i, xmin + i * dx);
		for (int i = 0; i <= M; i++)
			gsl_vector_set(tMesh, i, i * dt);
		wMatrix = gsl_matrix_alloc(N + 1, M + 1);
		//Apply Boundary Condition

		for (int j = 0; j <= N; j++)
			gsl_matrix_set(wMatrix, j, 0,
					GSL_MAX(exp(gsl_vector_get(xMesh, j)) - K, 0));
		for (int i = 0; i <= M; i++)
			gsl_matrix_set(wMatrix, 0, i, 0);
		for (int i = 0; i <= M; i++)
			gsl_matrix_set(wMatrix, N, i, 0);

		//calculate g vectors
		g1klamdan = vector<double>(N + 1);
		g1klamdap = vector<double>(N + 1);
		g2klamdan = vector<double>(N + 1);
		g2klamdap = vector<double>(N + 1);
		g2klamdan_plus1 = vector<double>(N + 1);
		g2klamdap_minus1 = vector<double>(N + 1);
		g2klamdap_plus1 = vector<double>(N + 1);

		//prestore g() value
		for (int k = 0; k <= N; k++) {
			g1klamdan[k] = g1(Y, k * dx * lamdan());
			g1klamdap[k] = g1(Y, k * dx * lamdap());
			g2klamdan[k] = g2(Y, k * dx * lamdan());
			g2klamdap[k] = g2(Y, k * dx * lamdap());
			g2klamdan_plus1[k] = g2(Y, k * dx * (lamdan() + 1));
			g2klamdap_minus1[k] = g2(Y, k * dx * (lamdap() - 1));
			g2klamdap_plus1[k] = g2(Y, k * dx * (lamdap() + 1));
		}

		//set A implicit three bands are diag,e,f (d is diagonal, e is upper diagonal, f is lower diagonal)
	}
	void solve() {

		//start the clock
		std::clock_t start;
		double duration;
		start = std::clock();

		gsl_vector *b = gsl_vector_alloc(N - 1);
		gsl_vector *diag = gsl_vector_alloc(N - 1);
		gsl_vector *e = gsl_vector_alloc(N - 2);
		gsl_vector *f = gsl_vector_alloc(N - 2);
		gsl_vector *x = gsl_vector_alloc(N - 1);
		double uvalue = u();
		double lvalue = l();
		for (int i = 0; i < N - 2; i++) {
			gsl_vector_set(e, i, uvalue);
			gsl_vector_set(f, i, lvalue);
			gsl_vector_set(diag, i, d(N, i + 1));
		}
		gsl_vector_set(diag, N - 2, d(N, N - 1));

		for (int j = 1; j <= M; j++) {
			//set  b
			for (int i = 1; i <= N - 1; i++) {
				double rhs;
				if (i == 1) {
					rhs = Wat(i, j - 1)
							+ dt * R(i, j - 1) / vega- lvalue*Wat(i-1,j-1);
					gsl_vector_set(b, i - 1, rhs);
				} else if (i == N - 1) {
					rhs = Wat(i, j - 1)
							+ dt * R(i, j - 1) / vega- uvalue*Wat(i+1,j-1);
					gsl_vector_set(b, i - 1, rhs);
				} else {
					rhs = Wat(i,j-1) + dt * R(i, j - 1) / vega;
					gsl_vector_set(b, i - 1, rhs);
				}
			}
			//solve AX = b
			gsl_linalg_solve_tridiag(diag, e, f, b, x);
			//update w matrix
			for (int i = 1; i <= N - 1; i++) {
				gsl_matrix_set(wMatrix, i, j, gsl_vector_get(x, i - 1));
			}
		}
		writematrix(wMatrix, N, M);
		cout << "M = " << M << " N=" << N << endl;
		cout << "UOC Premium:";
		cout << Wat(floor(N * ((log(S) - xmin) / (xmax - xmin))), M) << endl;
		duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
		cout << "Runtime(sec): " << duration << endl;

	}

};

int main() {
	int m = 0; int n= 0;
	cout<<"Please input M and N:"<<endl;
	cin>>m;
	cin>>n;
	PIDEsolver UOC(m,n);
	UOC.solve();

	return 0;
}
