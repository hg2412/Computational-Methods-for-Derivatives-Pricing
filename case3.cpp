/*
 * case3.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: Haoxiang
 */

#include <iostream>
#include <array>
#include <cmath>
#include <utility>
#include <vector>
#include <sstream>
#include <complex>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <fstream>
#include "simplex.h"

using namespace std;
typedef complex<double> dcomplex;

class Option {
public:
	string name;
	enum CP {
		call, put
	};
	double bid;
	double ask;
	double mid;
	double modelprice;
	double r;
	double q;
	double Spot;
	double maturity;
	CP cp;
	double K;

	void print() {
		cout << name << endl;
		cout << "Call or Put:" << cp << endl;
		cout << "bid:" << bid << endl;
		cout << "ask" << ask << endl;
		cout << "mid" << mid << endl;
		cout << "r" << r << endl;
		cout << "q" << q << endl;
		cout << "Spot" << Spot << endl;
		cout << "Maturity" << maturity << endl;
		cout << "Strike" << K << endl;
	}

	complex<double> HestonCF(dcomplex u, double t, double kappa, double theta,
			double sigma, double rho, double v0);
	complex<double> VGSACF(dcomplex u, double t, double sigma, double vega, double theta,
			double kappa, double eta, double lamda);
	double getfrFFTpriceHeston(double lamda, int n, double kappa, double theta,
			double sigma, double rho, double v0);
	double getfrFFTpriceVGSA(double lamda, int n, double sigma, double vega, double theta, double kappa, double eta2,double lamda2);
	double getBSprice(double sigma);
};

vector<Option> options;

inline complex<double> Option::HestonCF(dcomplex u, double t, double kappa,
		double theta, double sigma, double rho, double v0) {
	dcomplex gamma = sqrt(
			pow(sigma, 2) * (pow(u, 2) + u * dcomplex(0, 1))
					+ pow(kappa - rho * sigma * u * dcomplex(0, 1), 2));
	dcomplex phi = exp(
			u * dcomplex(0, 1) * log(Spot) + (r - q) * t * u * dcomplex(0, 1)
					+ kappa * theta * t
							* (kappa - rho * sigma * u * dcomplex(0, 1))
							/ pow(sigma, 2));
	phi = phi
			/ (pow(
					cosh(gamma * t * 0.5)
							+ (kappa - rho * sigma * u * dcomplex(0, 1)) / gamma
									* sinh(gamma * t * 0.5),
					2.0 * kappa * theta / pow(sigma, 2)));
	phi = phi
			* exp(
					-((pow(u, 2) + dcomplex(0, 1) * u) * v0)
							/ (gamma / tanh(gamma * t * 0.5) + kappa
									- rho * sigma * u * dcomplex(0, 1)));
	return phi;
}

dcomplex PhiVG(dcomplex u, double vega, double sigma, double theta){
	dcomplex t = 1.0- dcomplex(0,1) * u*theta * vega + sigma*sigma *u * u * vega/2.0;
	return -(1.0/vega) * log(t);
}

dcomplex Phi(dcomplex u, double t, double y0, double kappa, double eta, double lamda){
	dcomplex gamma = sqrt(kappa*kappa - 2* lamda*lamda * dcomplex(0,1)*u);
	dcomplex A = exp((kappa*kappa*eta*t)/(lamda*lamda)) / pow( cosh(gamma*t/2.0) + kappa/gamma*sinh(gamma*t/2.0), (2.0*kappa*eta)/(lamda*lamda));
	dcomplex B = 2.0*u*dcomplex(0,1)/(kappa + gamma/tanh(gamma*t/2.0));
	return (A * exp(B*y0));
}


inline complex<double >Option::VGSACF(dcomplex u, double t, double sigma, double vega, double theta,
		double kappa, double eta, double lamda){
	dcomplex t1 = exp(dcomplex(0,1)*u*(log(Spot) + (r-q)*t));
	dcomplex t2 = Phi(-dcomplex(0,1)*PhiVG(u,vega, sigma,theta),t,1/vega,kappa,eta,lamda);
	dcomplex t3 = pow(Phi(-dcomplex(0,1)*PhiVG(dcomplex(0,-1),vega, sigma,theta),t,1/vega,kappa,eta,lamda), u*dcomplex(0,1));

	return t1*t2/t3;
}

inline double Option::getBSprice(double sigma) {
	double S = Spot;
	double T = maturity / 360;
	double d1 = (log(S / K) + (r - q + sigma * sigma / 2.) * (T))
			/ (sigma * sqrt(T));
	double d2 = (log(S / K) + (r - q - sigma * sigma / 2.) * (T))
			/ (sigma * sqrt(T));
	double price = 0;
	if (cp == call)
		price = gsl_cdf_gaussian_P(d1, 1) * S * exp(-q * T)
				- gsl_cdf_gaussian_P(d2, 1) * K * exp(-r * T);
	if (cp == put)
		price = gsl_cdf_gaussian_P(-d2, 1) * K * exp(-r * T)
				- gsl_cdf_gaussian_P(-d1, 1) * S * exp(-q * T);
	return price;
}

inline double Option::getfrFFTpriceHeston(double lamda, int n, double kappa,
		double theta, double sigma, double rho, double v0) {
	int N = pow(2, n);
	array<double, 65536> x;
	array<double, 65536> Y;
	array<double, 65536> Z;
	array<double, 65537> C;
	array<double, 65536> KSI;
	double v = 0;
	double eta = 0.25;
	double alpha = 1;
	if (cp == put)
		alpha = -3;
	double gamma = eta * lamda / (2 * M_PI);
	double beta = log(K) - 0.5 * lamda * N;
	double T = maturity / 360;
	gsl_complex z = gsl_complex_rect(0, 0);
	dcomplex zdcomp(0, 0);
	//generate X
	for (int i = 1; i <= N; i++) {
		v = (i - 1) * eta;
		zdcomp = HestonCF(dcomplex(v, -(1 + alpha)), T, kappa, theta, sigma,
				rho, v0);
		z = gsl_complex_rect(zdcomp.real(), zdcomp.imag());
		z = gsl_complex_mul(z,
				gsl_complex_div(gsl_complex_rect(eta * exp(-r * T), 0),
						gsl_complex_mul(gsl_complex_rect(alpha, v),
								gsl_complex_rect(alpha + 1, v))));
		z = gsl_complex_mul(z, gsl_complex_exp(gsl_complex_rect(0, -beta * v)));
		if (i == 1) {
			x.at(2 * i - 2) = z.dat[0] * 0.5;
			x.at(2 * i - 1) = z.dat[1] * 0.5;
		} else {
			x.at(2 * i - 2) = z.dat[0];
			x.at(2 * i - 1) = z.dat[1];
		}
	}

	//set vector Y and Z
	gsl_complex xcomp;
	gsl_complex ycomp;
	gsl_complex zcomp;
	for (int i = 1; i <= N; i++) {
		xcomp = gsl_complex_rect(x.at(2 * i - 2), x.at(2 * i - 1));
		ycomp = gsl_complex_mul(
				gsl_complex_exp(
						gsl_complex_rect(0, -M_PI * gamma * (i - 1) * (i - 1))),
				xcomp);
		zcomp = gsl_complex_exp(
				gsl_complex_rect(0, M_PI * gamma * (i - 1) * (i - 1)));
		Y.at(2 * i - 2) = ycomp.dat[0];
		Y.at(2 * i - 1) = ycomp.dat[1];
		Z.at(2 * i - 2) = zcomp.dat[0];
		Z.at(2 * i - 1) = zcomp.dat[1];
		Z.at(2 * (2 * N - i + 1) - 2) = zcomp.dat[0];
		Z.at(2 * (2 * N - i + 1) - 1) = zcomp.dat[1];
	}
	for (int i = N + 1; i <= 2 * N; i++) {
		Y.at(2 * i - 2) = 0;
		Y.at(2 * i - 1) = 0;
	}
	//FFT on Y and Z
	gsl_fft_complex_radix2_forward(Y.data(), 1, 2 * N);
	gsl_fft_complex_radix2_forward(Z.data(), 1, 2 * N);
	gsl_complex ksi;
	for (int i = 1; i <= 2 * N; i++) {
		ksi = gsl_complex_mul(
				gsl_complex_rect(Y.at(2 * i - 2), Y.at(2 * i - 1)),
				gsl_complex_rect(Z.at(2 * i - 2), Z.at(2 * i - 1)));
		KSI.at(2 * i - 2) = ksi.dat[0];
		KSI.at(2 * i - 1) = ksi.dat[1];
	}
	//inverse FFT on Ksi vector
	gsl_fft_complex_radix2_inverse(KSI.data(), 1, 2 * N);
	//calculate Call prices
	for (int i = 1; i <= N; i++) {
		ksi = gsl_complex_rect(KSI.at(2 * i - 2), KSI.at(2 * i - 1));
		C.at(i) = exp(-alpha * (log(K) - (N / 2 - i + 1) * lamda)) / M_PI;
		C.at(i) =
				C.at(
						i) *GSL_REAL( gsl_complex_mul(ksi, gsl_complex_exp(gsl_complex_rect(0,-M_PI*gamma*(i-1)*(i-1)))));
	}

	//output option price at given Strike price
	double *cp = C.data();
	return cp[N / 2 + 1];
}


inline double Option::getfrFFTpriceVGSA(double lamda, int n,  double sigma, double vega, double theta, double kappa, double eta2,double lamda2) {
	int N = pow(2, n);
	array<double, 65536> x;
	array<double, 65536> Y;
	array<double, 65536> Z;
	array<double, 65537> C;
	array<double, 65536> KSI;
	double v = 0;
	double eta = 0.25;
	double alpha = 1;
	if (cp == put)
		alpha = -3;
	double gamma = eta * lamda / (2 * M_PI);
	double beta = log(K) - 0.5 * lamda * N;
	double T = maturity / 360;
	gsl_complex z = gsl_complex_rect(0, 0);
	dcomplex zdcomp(0, 0);
	//generate X
	for (int i = 1; i <= N; i++) {
		v = (i - 1) * eta;
		zdcomp = VGSACF(dcomplex(v, -(1 + alpha)), T, sigma, vega, theta,kappa,eta2, lamda2);
		z = gsl_complex_rect(zdcomp.real(), zdcomp.imag());
		z = gsl_complex_mul(z,
				gsl_complex_div(gsl_complex_rect(eta * exp(-r * T), 0),
						gsl_complex_mul(gsl_complex_rect(alpha, v),
								gsl_complex_rect(alpha + 1, v))));
		z = gsl_complex_mul(z, gsl_complex_exp(gsl_complex_rect(0, -beta * v)));
		if (i == 1) {
			x.at(2 * i - 2) = z.dat[0] * 0.5;
			x.at(2 * i - 1) = z.dat[1] * 0.5;
		} else {
			x.at(2 * i - 2) = z.dat[0];
			x.at(2 * i - 1) = z.dat[1];
		}
	}

	//set vector Y and Z
	gsl_complex xcomp;
	gsl_complex ycomp;
	gsl_complex zcomp;
	for (int i = 1; i <= N; i++) {
		xcomp = gsl_complex_rect(x.at(2 * i - 2), x.at(2 * i - 1));
		ycomp = gsl_complex_mul(
				gsl_complex_exp(
						gsl_complex_rect(0, -M_PI * gamma * (i - 1) * (i - 1))),
				xcomp);
		zcomp = gsl_complex_exp(
				gsl_complex_rect(0, M_PI * gamma * (i - 1) * (i - 1)));
		Y.at(2 * i - 2) = ycomp.dat[0];
		Y.at(2 * i - 1) = ycomp.dat[1];
		Z.at(2 * i - 2) = zcomp.dat[0];
		Z.at(2 * i - 1) = zcomp.dat[1];
		Z.at(2 * (2 * N - i + 1) - 2) = zcomp.dat[0];
		Z.at(2 * (2 * N - i + 1) - 1) = zcomp.dat[1];
	}
	for (int i = N + 1; i <= 2 * N; i++) {
		Y.at(2 * i - 2) = 0;
		Y.at(2 * i - 1) = 0;
	}
	//FFT on Y and Z
	gsl_fft_complex_radix2_forward(Y.data(), 1, 2 * N);
	gsl_fft_complex_radix2_forward(Z.data(), 1, 2 * N);
	gsl_complex ksi;
	for (int i = 1; i <= 2 * N; i++) {
		ksi = gsl_complex_mul(
				gsl_complex_rect(Y.at(2 * i - 2), Y.at(2 * i - 1)),
				gsl_complex_rect(Z.at(2 * i - 2), Z.at(2 * i - 1)));
		KSI.at(2 * i - 2) = ksi.dat[0];
		KSI.at(2 * i - 1) = ksi.dat[1];
	}
	//inverse FFT on Ksi vector
	gsl_fft_complex_radix2_inverse(KSI.data(), 1, 2 * N);
	//calculate Call prices
	for (int i = 1; i <= N; i++) {
		ksi = gsl_complex_rect(KSI.at(2 * i - 2), KSI.at(2 * i - 1));
		C.at(i) = exp(-alpha * (log(K) - (N / 2 - i + 1) * lamda)) / M_PI;
		C.at(i) =
				C.at(
						i) *GSL_REAL( gsl_complex_mul(ksi, gsl_complex_exp(gsl_complex_rect(0,-M_PI*gamma*(i-1)*(i-1)))));
	}

	//output option price at given Strike price
	double *cp = C.data();
	return cp[N / 2 + 1];
}




//relative error squraed weighted inversely by bid-ask spread
double HestonLossFunction1(vector<Option> &options, double kappa, double theta,
		double sigma, double rho, double v0) {
	double SSE = 0;
	for (auto &opt : options) {
		SSE += pow(
				(opt.mid
						- opt.getfrFFTpriceHeston(0.1, 8, kappa, theta, sigma, rho,
								v0)), 2);
	}
	return SSE;

}

//relative error squared weighted inversely by bid-ask spread
double HestonLossFunction2(vector<Option> &options, double kappa, double theta,
		double sigma, double rho, double v0) {
	double SSE = 0;
	for (auto &opt : options) {
		SSE += pow(
				(opt.mid
						- opt.getfrFFTpriceHeston(0.1, 8, kappa, theta, sigma, rho,
								v0)), 2)
				/ (opt.ask - opt.bid);
	}
	return SSE;
}

template<class Con>
void printcon(const Con& c) {
	std::cout.precision(8);
	cout << "results are:";
	copy(c.begin(), c.end(),
			ostream_iterator<typename Con::value_type>(cout, ","));
	cout << endl;
}

double Heston(vector<double> x) {
	return HestonLossFunction1(options, x[0], x[1], x[2], x[3], x[4]);
}

double Heston2(vector<double> x) {
	return HestonLossFunction2(options, x[0], x[1], x[2], x[3], x[4]);
}

//A function read csv to a vector of options
void readCSV(vector<Option> &options, string filename) {
	string line, csvItem;
	ifstream myfile(filename);
	int lineNumber = 0;
	if (myfile.is_open()) {
		while (getline(myfile, line, '\r')) {
			lineNumber++;
			if (lineNumber >= 2) {
				istringstream myline(line);
				int colnum = 0;
				Option opt;
				while (getline(myline, csvItem, ',')) {
					colnum++;
					switch (colnum) {
					case 1:
						opt.name = csvItem;
						break;
					case 2:
						opt.bid = stod(csvItem);
						break;
					case 3:
						opt.ask = stod(csvItem);
						break;
					case 4:
						opt.r = stod(csvItem) * 0.01;
						break;
					case 5:
						opt.q = stod(csvItem) * 0.01;
						break;
					case 7:
						opt.Spot = stod(csvItem);
						break;
					case 9:
						opt.cp = (csvItem == "C" ? Option::call : Option::put);
						break;
					case 10:
						opt.K = stod(csvItem);
						break;
					case 11:
						opt.maturity = stod(csvItem);
						break;
					default:
						;
					}
				}
				opt.mid = 0.5 * (opt.bid + opt.ask);
				options.push_back(opt);

			}
		}
		myfile.close();
	}
	return;

}

void CalibrationHeston1 (vector<Option> &options){
	vector<double> init;
	init.push_back(5);
	init.push_back(0.1);
	init.push_back(1);
	init.push_back(-0.6);
	init.push_back(0.06);

	double a0[] = { 4.5, 0.1, 0.5, -0.3, 0.01 };
	double a1[] = { 4.3, 0.2, 0.7, -0.5, 0.1 };
	double a2[] = { 4.8, 0.15, 1.1, -0.1, 0.2 };
	double a3[] = { 5.1, 0.05, 1.3, 0.5, 0.25 };
	double a4[] = { 5.2, 0.17, 1  , 0.3, 0.05 };
	double a5[] = { 5.5, 0.07, 0.9, 0.1, 0.15 };

	vector<vector<double> > simplex;
	simplex.push_back(vector<double>(a0, a0 + 5));
	simplex.push_back(vector<double>(a1, a1 + 5));
	simplex.push_back(vector<double>(a2, a2 + 5));
	simplex.push_back(vector<double>(a3, a3 + 5));
	simplex.push_back(vector<double>(a4, a4 + 5));
	simplex.push_back(vector<double>(a5, a5 + 5));

	using BT::Simplex;
	cout << Heston(init);

	//result: -0.0030718125  -0.17335788  0.61328266  -0.11075356  -0.056074466
	//result:  0.49376  -0.0702144  0.44592  -0.446912  -0.0149824
//	cout << "Function achieves minimum at:" << endl;
	auto result = Simplex(Heston, init,1e-3, simplex);
	cout<<"Equally Weighted Result:"<<endl;
	cout<<"Squared Relative Error"<<Heston(result)<<endl;
	cout<<"Optimal Parameter:";
	printcon(result);

	ofstream fout("Calibration Heston 1.txt");
	fout<<"Squared Error:"<<endl;
	fout<<Heston(result)<<endl;
	fout<<"Optimal Parameters:";
	for (auto &r:result)
		fout<<r<<" ";
	fout.close();
}


void CalibrationHeston2 (vector<Option> &options){
	vector<double> init;
	init.push_back(5);
	init.push_back(0.1);
	init.push_back(1);
	init.push_back(-0.6);
	init.push_back(0.06);

	double a0[] = { 4.5, 0.1, 0.5, -0.7, 0.01 };
	double a1[] = { 4.3, 0.2, 0.7, -0.5, 0.1 };
	double a2[] = { 4.8, 0.15, 1.1, -0.1, 0.2 };
	double a3[] = { 5.1, 0.05, 1.3, -0.2, 0.25 };
	double a4[] = { 5.2, 0.17, 1 , 0.2, 0.05 };
	double a5[] = { 5.5, 0.07, 0.9, 0.1, 0.15 };

	vector<vector<double> > simplex;
	simplex.push_back(vector<double>(a0, a0 + 5));
	simplex.push_back(vector<double>(a1, a1 + 5));
	simplex.push_back(vector<double>(a2, a2 + 5));
	simplex.push_back(vector<double>(a3, a3 + 5));
	simplex.push_back(vector<double>(a4, a4 + 5));
	simplex.push_back(vector<double>(a5, a5 + 5));

	using BT::Simplex;
	cout << Heston2(init);


	auto result = Simplex(Heston2, init, 1e-3, simplex);

	cout<<"Weighted inversely by spread"<<endl;
	cout<<"Squared Relative Error"<<Heston2(result)<<endl;
	cout<<"Optimal Parameter:";
	printcon(result);


	ofstream fout("Calibration Heston 2.txt");
	fout<<"Squared Error:"<<endl;
	fout<<Heston2(result)<<endl;
	fout<<"Optimal Parameters:";
	for (auto &r:result)
		fout<<r<<" ";
	fout.close();

}

//relative error squared equally weighted
double VGSALossFunction1(vector<Option> &options,double sigma, double vega, double theta, double kappa, double eta2,double lamda2 ) {
	double SSE = 0;
	for (auto &opt : options) {
		SSE += pow(
				(opt.mid - opt.getfrFFTpriceVGSA(0.1, 8, sigma, vega, theta,kappa,eta2,lamda2)), 2);
	}
	return SSE;

}

//relative error squared weighted inversely by bid-ask spread
double VGSALossFunction2(vector<Option> &options, double sigma, double vega, double theta, double kappa, double eta2,double lamda2) {
	double SSE = 0;
	for (auto &opt : options) {
		SSE += pow(
				(opt.mid
						- opt.getfrFFTpriceVGSA(0.1, 8, sigma, vega, theta, kappa, eta2, lamda2)), 2)
				/ (opt.ask - opt.bid);
	}
	return SSE;
}

double VGSA1(vector<double> x) {
	return VGSALossFunction1(options, x[0], x[1], x[2], x[3], x[4],x[5]);
}

double VGSA2(vector<double> x) {
	return VGSALossFunction2(options, x[0], x[1], x[2], x[3], x[4],x[5]);
}

void CalibrationVGSA1 (vector<Option> &options){
	vector<double> init;
	init.push_back(0.1022);
	init.push_back(0.1819);
	init.push_back(-0.0761);
	init.push_back(8.1143);
	init.push_back(2.8060);
	init.push_back(10.36);

	double a0[] = { 0.05, 0.1, -0.1, 7, 3 ,10 };
	double a1[] = { 0.1, 0.2, -0.2, 8, 2.5 , 11 };
	double a2[] = { 0.08, 0.3, 0.2, 8.5, 4 ,9.5 };
	double a3[] = { 0.06, 0.22, -0.05,9, 3.5,  11.5 };
	double a4[] = { 0.12, 0.15, -0.07,6.5, 2.7 , 10.5 };
	double a5[] = { 0.13, 0.25, -0.08, 9.5, 3.1, 9 };
	double a6[] = { 0.14, 0.3, 0.1, 8.7, 3.4, 9.8 };

	vector<vector<double> > simplex;
	simplex.push_back(vector<double>(a0, a0 + 6));
	simplex.push_back(vector<double>(a1, a1 + 6));
	simplex.push_back(vector<double>(a2, a2 + 6));
	simplex.push_back(vector<double>(a3, a3 + 6));
	simplex.push_back(vector<double>(a4, a4 + 6));
	simplex.push_back(vector<double>(a5, a5 + 6));
	simplex.push_back(vector<double>(a6, a6 + 6));

	using BT::Simplex;
	cout << VGSA1(init);

	//result: -0.0030718125  -0.17335788  0.61328266  -0.11075356  -0.056074466
	//result:  0.49376  -0.0702144  0.44592  -0.446912  -0.0149824
//	cout << "Function achieves minimum at:" << endl;
	auto result = Simplex(VGSA1, init, 1e-3, simplex);
	cout<<"Equally Weighted Result:"<<endl;
	cout<<"Squared Relative Error"<<VGSA1(result)<<endl;
	cout<<"Optimal Parameter:";
	printcon(result);

	ofstream fout("Calibration VGSA 1.txt");
	fout<<"Squared Error:"<<endl;
	fout<<VGSA1(result)<<endl;
	fout<<"Optimal Parameters:";
	for (auto &r:result)
		fout<<r<<" ";
	fout.close();
}

void CalibrationVGSA2 (vector<Option> &options){
	vector<double> init;
	init.push_back(0.1022);
	init.push_back(0.1819);
	init.push_back(-0.0761);
	init.push_back(8.1143);
	init.push_back(2.8060);
	init.push_back(10.3646);

	double a0[] = { 0.05, 0.1, -0.1, 7, 5 , 10};
	double a1[] = { 0.1, 0.2, -0.25, 5, 8, 3, 8 };
	double a2[] = { 0.08, 0.3, 0.2, 2, 1, 13 };
	double a3[] = { 0.06, 0.4, -0.5,6, 3 , 5};
	double a4[] = { 0.12, 0.15, -0.3,8, 7 ,7};
	double a5[] = { 0.13, 0.5, -0.15, 4, 6, 9 };
	double a6[] = { 0.14, 0.33, 0.1, 1, 2,13 };

	vector<vector<double> > simplex;
	simplex.push_back(vector<double>(a0, a0 + 6));
	simplex.push_back(vector<double>(a1, a1 + 6));
	simplex.push_back(vector<double>(a2, a2 + 6));
	simplex.push_back(vector<double>(a3, a3 + 6));
	simplex.push_back(vector<double>(a4, a4 + 6));
	simplex.push_back(vector<double>(a5, a5 + 6));
	simplex.push_back(vector<double>(a6, a6 + 6));

	using BT::Simplex;
	cout << VGSA2(init);

	//result: -0.0030718125  -0.17335788  0.61328266  -0.11075356  -0.056074466
	//result:  0.49376  -0.0702144  0.44592  -0.446912  -0.0149824
//	cout << "Function achieves minimum at:" << endl;
	auto result = Simplex(VGSA2, init, 1e-3, simplex);
	cout<<"Inversely Weighted Result:"<<endl;
	cout<<"Squared Relative Error"<<VGSA2(result)<<endl;
	cout<<"Optimal Parameter:";
	printcon(result);

	ofstream fout("Calibration VGSA 2.txt");
	fout<<"Squared Error:"<<endl;
	fout<<VGSA2(result)<<endl;
	fout<<"Optimal Parameters:";
	for (auto &r:result)
		fout<<r<<" ";
	fout.close();
}

void HestonSurface(double spot,double r, double q, int N, int M, double kappa, double theta,
		double sigma, double rho, double v0 ){
	double dt = 1;
	double dk = 1;
	vector<double> T;
	vector<double> K;
	vector< vector<double> > opt_matrix(M,vector<double>(N));
	vector< vector<double> > sigma_matrix(M,vector<double>(N));
	for (int i = 0;i <M;i++){
		T.push_back(i);
	}
	for (int j = 0;j <N;j++){
		K.push_back(spot - N/2 + j);
	}

	for (int i=0;i<M;i++){
		for (int j=0; j<N;j++){
			Option opt;
			opt.K = K[j];
			opt.maturity =T[i];
			opt.r = r;
			opt.q = q;
			opt.cp = Option::call;
			opt.Spot = spot;
			opt_matrix[i][j] = opt.getfrFFTpriceHeston(0.1,8, kappa,theta,sigma,rho,v0);
		}
	}

	for (int i= 1;i<M-1;i++){
		for (int j=1;j<N-1;j++){
			double C = opt_matrix[i][j];
			double k = K[j];
			double dCT = 180*(opt_matrix[i+1][j] - opt_matrix[i-1][j]);
			double dCK = 0.5*(opt_matrix[i][j+1] - opt_matrix[i][j-1]);
			double ddCK =  (opt_matrix[i][j+1] + opt_matrix[i][j-1] - 2* opt_matrix[i][j]);
			sigma_matrix[i][j] =sqrt(2*(dCT+q*C + (r-q)*k*dCK)/(k*k*ddCK));
		}
	}


	ofstream fout("Heston Vol Surface.csv");
	for (int i = 1; i < M-1; i++) {
		for (int j = 1; j < N-1; j++)
			fout << sigma_matrix[i][j] << ",";
		fout << "\n";
	}

	ofstream f1("Heston T.csv");
	ofstream f2("Heston K.csv");

	for (int i = 1; i < M-1; i++) {
		f1 <<T[i]<<",";
	}

	for (int i = 1; i < N-1; i++) {
			f2 <<K[i]<<",";
	}
	fout.close();
	f1.close();
	f2.close();

}


void VGSASurface(double spot,double r, double q, int N, int M, double sigma, double vega, double theta, double kappa, double eta2,double lamda2 ){
	double dt = 1;
	double dk = 1;
	vector<double> T;
	vector<double> K;
	vector< vector<double> > opt_matrix(M,vector<double>(N));
	vector< vector<double> > sigma_matrix(M,vector<double>(N));
	for (int i = 0;i <M;i++){
		T.push_back(i);
	}
	for (int j = 0;j <N;j++){
		K.push_back(spot - N/2 + j);
	}

	for (int i=0;i<M;i++){
		for (int j=0; j<N;j++){
			Option opt;
			opt.K = K[j];
			opt.maturity =T[i];
			opt.r = r;
			opt.q = q;
			opt.cp = Option::call;
			opt.Spot = spot;
			opt_matrix[i][j] = opt.getfrFFTpriceVGSA(0.1, 8, sigma,vega,theta,kappa,eta2,lamda2);
		}
	}

	for (int i= 1;i<M-1;i++){
		for (int j=1;j<N-1;j++){
			double C = opt_matrix[i][j];
			double k = K[j];
			double dCT = 180*(opt_matrix[i+1][j] - opt_matrix[i-1][j]);
			double dCK = 0.5*(opt_matrix[i][j+1] - opt_matrix[i][j-1]);
			double ddCK =  (opt_matrix[i][j+1] + opt_matrix[i][j-1] - 2* opt_matrix[i][j]);
			sigma_matrix[i][j] =sqrt(2*(dCT+q*C + (r-q)*k*dCK)/(k*k*ddCK));
		}
	}


	ofstream fout("VGSA Vol Surface.csv");

	ofstream f1("VGSA T.csv");
	ofstream f2("VGSA K.csv");

	for (int i = 1; i < M-1; i++) {
		for (int j = 1; j < N-1; j++)
			fout << sigma_matrix[i][j] << ",";
		fout << "\n";
	}

	for (int i = 1; i < M-1; i++) {
		f1 <<T[i]<<",";
	}

	for (int i = 1; i < N-1; i++) {
			f2 <<K[i]<<",";
	}

	fout.close();
	f1.close();
	f2.close();

}



int main() {
	readCSV(options, "sheet4.csv");
/*	Option opt;
	opt.Spot = 1900;
	opt.maturity = 90;
	opt.r = 0.0025;
	opt.q = 0.0187;
	opt.K = 2100;
	opt.cp = Option::call;
	cout << "Heston Price:"
			<< opt.getfrFFTpriceHeston(0.25, 8, 4, 0.06, 0.5, -0.5, 0.06) << endl;

	cout << "VGSA Price:"
			<< opt.getfrFFTpriceVGSA(0.01, 8, 0.1022,0.1819,-0.0761,8.1143,2.8060,10.3046) << endl;

	cout << "BS price" << opt.getBSprice(0.06) << endl;
*/
//	CalibrationHeston1(options); // result: 3.8,0.063,1.5,0.04,0.023,
//	CalibrationHeston2(options); // result: 3.8,0.063,1.5,0.04,0.023,
//	CalibrationVGSA1(options);//-0.0002315724,0.030320143,-0.033711403,14.275396,6.583543,22.244003,
//	CalibrationVGSA2(options); // -8.5145141e-07,1.4347281e-06,0.00017082866,9.3250387,-3.4004209,19.968406,
	HestonSurface(options[0].Spot, options[0].r, options[0].q, 202,202, 4.47401,0.0359928,0.46629,-0.92927,0.0143429 );
	VGSASurface(options[0].Spot, options[0].r, options[0].q, 202,202,0.024863,0.065975,-0.0517436,8.34739,3.98845,10.3674);
}

