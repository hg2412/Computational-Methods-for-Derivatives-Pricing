#include <iostream>
#include <array>
#include <cmath>
#include <utility>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

using namespace std;

void print_complex(gsl_complex z){
	cout<<"Real:"<<z.dat[0]<<" Imag:"<<z.dat[1]<<endl;
}

class PricingEngine{
public:
	double * getFFTprice();
	double * getfrFFTprice(double);
	double getBSprice();
	double getCosPricing(double a, double b);
	pair<double, double> getCosRange(double L);
	double cos_phi(double a, double b, double c, double d,int k);
	double cos_chi(double a, double b, double c, double d, int k);
	double cos_A(double a, double b, int k);
	gsl_complex cfunction(gsl_complex u); //characteristic function of log(S)
	gsl_complex cfunction2(gsl_complex u); //characteristic function of log(S/K)
	PricingEngine(){
		K = 2000; //K = 2000, 2100, 2200
		S = 1900;
		T = 0.25;
		sigma = 0.36;
		r = 0.0025;
		q = 0.0187;
		eta = 0.25;//warning! changed from 0.25
		alpha = 1;
		n = 8;//warning! has been changed from 9
	}

	double getAlpha() const {
		return alpha;
	}

	void setAlpha(double alpha) {
		this->alpha = alpha;
	}

	double getK() const {
		return K;
	}

	void setK(double k) {
		K = k;
	}

	int getN() const {
		return n;
	}

	void setN(int n) {
		this->n = n;
	}

private:
	double K;  //Strike Price
	double S;  //Stock Price
	double T;  //Time to Maturity
	double sigma; //volatility
	double r;  //risk-free rate
	double q;  //dividend rate
	int n;   //N = 2^n
	double alpha;  //damping factor
	double eta;
};

inline double * PricingEngine::getFFTprice() {
	int N = pow(2,n);
	double lamda = (2 * M_PI) / (N * eta);
	array<double, 65536> x;
	double v = 0;
	double beta = log(K)-0.5*lamda*N;
	array<double, 65536> C;

	gsl_complex z = gsl_complex_rect(0,0);
	//generate X
	for (int i=1;i<=N;i++){
		v = (i-1)*eta;
		z = cfunction(gsl_complex_rect(v,-(1+alpha)));
		z = gsl_complex_mul(z, gsl_complex_div(gsl_complex_rect(eta*exp(-r*T),0), gsl_complex_mul( gsl_complex_rect(alpha,v) ,gsl_complex_rect(alpha+1,v))));
		z = gsl_complex_mul(z, gsl_complex_exp(gsl_complex_rect(0,-beta*v)));
		if (i == 1){
			x.at(2*i-2) = z.dat[0]*0.5;
			x.at(2*i-1) = z.dat[1]*0.5;
		}else{
			x.at(2*i-2) = z.dat[0];
			x.at(2*i-1) = z.dat[1];
		}
	}
	//FFT
	gsl_fft_complex_radix2_forward(x.data(), 1,N );
	//get Call price
	for (int i=1;i<=N;i++){
		C.at(i) = x.at(2*i-2)*exp(-alpha*(log(K)-(N/2 - i+1)*lamda))/M_PI;
	}
	double *cp = C.data();
	cout<<cp[N/2+1]<<" ";
	return C.data();
}
//frFFT pricing
inline double* PricingEngine::getfrFFTprice(double lamda) {
	int N = pow(2,n);
	array<double, 65536> x;
	array<double, 65536> Y;
	array<double, 65536> Z;
	array<double, 65537> C;
	array<double, 65536> KSI;
	double v = 0;
	double gamma = eta*lamda/(2*M_PI);
	double beta = log(K)-0.5*lamda*N;
	gsl_complex z = gsl_complex_rect(0,0);
	//generate X
	for (int i=1;i<=N;i++){
		v = (i-1)*eta;
		z = cfunction(gsl_complex_rect(v,-(1+alpha)));
		z = gsl_complex_mul(z, gsl_complex_div(gsl_complex_rect(eta*exp(-r*T),0), gsl_complex_mul( gsl_complex_rect(alpha,v) ,gsl_complex_rect(alpha+1,v))));
		z = gsl_complex_mul(z, gsl_complex_exp(gsl_complex_rect(0,-beta*v)));
		if (i == 1){
			x.at(2*i-2) = z.dat[0]*0.5;
			x.at(2*i-1) = z.dat[1]*0.5;
		}else{
			x.at(2*i-2) = z.dat[0];
			x.at(2*i-1) = z.dat[1];
		}
	}

	//set vector Y and Z
	gsl_complex xcomp;
	gsl_complex ycomp;
	gsl_complex zcomp;
	for (int i =1;i<=N;i++){
		xcomp = gsl_complex_rect(x.at(2*i-2),x.at(2*i-1));
		ycomp = gsl_complex_mul(gsl_complex_exp(gsl_complex_rect(0,-M_PI*gamma*(i-1)*(i-1))),xcomp);
		zcomp = gsl_complex_exp(gsl_complex_rect(0,M_PI*gamma*(i-1)*(i-1)));
		Y.at(2*i-2) = ycomp.dat[0];
		Y.at(2*i-1) = ycomp.dat[1];
		Z.at(2*i-2) = zcomp.dat[0];
		Z.at(2*i-1) = zcomp.dat[1];
		Z.at(2*(2*N-i+1)-2) = zcomp.dat[0];
		Z.at(2*(2*N-i+1)-1) = zcomp.dat[1];
	}
	for (int i=N+1;i<=2*N;i++){
		Y.at(2*i-2) = 0;
		Y.at(2*i-1) = 0;
	}
	//FFT on Y and Z
	gsl_fft_complex_radix2_forward(Y.data(), 1,2*N );
	gsl_fft_complex_radix2_forward(Z.data(), 1,2*N );
	gsl_complex ksi;
	for (int i =1;i<=2*N;i++){
		ksi = gsl_complex_mul(gsl_complex_rect(Y.at(2*i-2),Y.at(2*i-1)),gsl_complex_rect(Z.at(2*i-2),Z.at(2*i-1)));
		KSI.at(2*i-2) = ksi.dat[0];
		KSI.at(2*i-1) = ksi.dat[1];
	}
	//inverse FFT on Ksi vector
	gsl_fft_complex_radix2_inverse(KSI.data(), 1,2*N );
	//calculate Call prices
	for (int i=1;i<=N;i++){
		ksi = gsl_complex_rect(KSI.at(2*i-2),KSI.at(2*i-1));
		C.at(i) = exp(-alpha*(log(K)-(N/2 - i+1)*lamda))/M_PI;
		C.at(i) = C.at(i) *GSL_REAL( gsl_complex_mul(ksi, gsl_complex_exp(gsl_complex_rect(0,-M_PI*gamma*(i-1)*(i-1)))));
	}

	//output option price at given Strike price
	double *cp = C.data();
	cout<<cp[N/2+1]<<" ";
	return C.data();
}
//BS formula pricing
inline double PricingEngine::getBSprice(){
	double d1=(log(S/K) + (r -q + sigma*sigma/2.)*(T))/(sigma*sqrt(T));
	double d2=(log(S/K) + (r -q - sigma*sigma/2.)*(T))/(sigma*sqrt(T));
	double price = gsl_cdf_gaussian_P(d1,1)*S*exp(-q*T) - gsl_cdf_gaussian_P(d2,1)*K*exp(-r*T);
	cout<<"BS Price:"<<endl;
	cout<<price<<endl;
	return price;
}
//phi function for cosine method
inline double PricingEngine::cos_phi(double a, double b, double c, double d, int k) {
	if (k == 0)
		return (d - c);
	double phi = sin(k*M_PI*(d-a)/(b-a)) - sin(k*M_PI*(c-a)/(b-a));
	phi = phi * (b-a)/(k*M_PI);
	return phi;
}
//chi function for cosine method
inline double PricingEngine::cos_chi(double a, double b, double c, double d, int k) {
	double chi1 = 1/ ( 1 + pow( ( k*M_PI / (b-a) ),2) );
	double chi2 = cos(k*M_PI*(d-a)/(b-a))*exp(d) - cos(k*M_PI*(c-a)/(b-a)*exp(c));
	chi2 = chi2 + (k*M_PI / (b-a)) * sin(k*M_PI*(d-a)/(b-a))*exp(d) - k*M_PI / (b-a) * sin(k*M_PI*(c-a)/(b-a))*exp(c);
	return chi1*chi2;
}
//Cosine method pricing
inline double PricingEngine::getCosPricing(double a, double b) {
	double V = 0;
	double vk = 0;
	int N = pow(2,n);
	for (int k=0;k<N;k++){
		vk = K*(2/(b-a))*(cos_chi(a,b,0,b,k) - cos_phi(a,b,0,b,k));
		if (k == 0) vk = vk * 0.5;
		V = V + vk*cos_A(a,b,k);
	}
	V = V*exp(-r*T)*(b-a)/2;
	cout<<V<<" ";
	return V;
}

//generate A coefficients of cosine method
inline double PricingEngine::cos_A(double a, double b, int k) {
	double term1 = k*M_PI/(b-a);
	gsl_complex comp;
	comp = gsl_complex_mul(cfunction2(gsl_complex_rect(term1,0)), gsl_complex_exp(gsl_complex_rect(0,-a*term1)));
	double A = GSL_REAL(comp)*2/(b-a);
	return A;
}
//characteristic function of logS
inline gsl_complex PricingEngine::cfunction(gsl_complex u) {
	gsl_complex a1 = gsl_complex_mul_imag(u , log(S) + (r - q - 0.5*sigma*sigma)*T);
	gsl_complex a2 = gsl_complex_mul_real(gsl_complex_mul(u,u) , - 0.5*sigma*sigma*T);
	return gsl_complex_exp(gsl_complex_add(a1,a2));
}
//generate the suggested range for Cosine method
inline pair<double, double> PricingEngine::getCosRange(double L) {
	double a = (r-q)*T - L*sqrt(pow(sigma,2)*T);
	double b = (r-q)*T + L*sqrt(pow(sigma,2)*T);
	return pair<double, double>(a,b);
}

inline gsl_complex PricingEngine::cfunction2(gsl_complex u) {
	gsl_complex a1 = gsl_complex_mul_imag(u , log(S/K) + (r - q - 0.5*sigma*sigma)*T);
	gsl_complex a2 = gsl_complex_mul_real(gsl_complex_mul(u,u) , - 0.5*sigma*sigma*T);
	return gsl_complex_exp(gsl_complex_add(a1,a2));
}

int main()
{
	PricingEngine OptionPricing;
	double K[3] = {2000,2100,2200};
	double Alpha[4] = {0.4,1.0,1.4,3.0};
	int ns[4] = {9,11,13,15};
	cout<<"FFT"<<endl;
	for (int i=0;i<3;i++){
		cout<<"K="<<K[i]<<endl;
		for (int j=0;j<4;j++){
			cout<<"Alpha="<<Alpha[j]<<endl;
			for (int k=0;k<4;k++){
				OptionPricing.setK(K[i]);
				OptionPricing.setAlpha(Alpha[j]);
				OptionPricing.setN(ns[k]);
				OptionPricing.getFFTprice();
			}
			cout<<endl;
		}
	}

	int nfrfft[4] = {6,7,8,9};
	cout<<"frFFT"<<endl;
	for (int i=0;i<3;i++){
		cout<<"K="<<K[i]<<endl;
		for (int j=0;j<4;j++){
			cout<<"Alpha="<<Alpha[j]<<endl;
			for (int k=0;k<4;k++){
				OptionPricing.setK(K[i]);
				OptionPricing.setAlpha(Alpha[j]);
				OptionPricing.setN(nfrfft[k]);
				OptionPricing.getfrFFTprice(0.1);
			}
			cout<<endl;
		}
	}
	cout<<"COS"<<endl;
	OptionPricing.setN(6);
	cout<<"n=6"<<endl;
	double bound[4] = {1,4,8,12};
	for (int i=0;i<3;i++){
		cout<<"K="<<K[i]<<endl;
		for (int j=0;j<4;j++){
			OptionPricing.setK(K[i]);
			OptionPricing.getCosPricing(-bound[j],bound[j]);
		}
		cout<<endl;
	}

	OptionPricing.setN(9);
	cout<<"n=9"<<endl;
	for (int i=0;i<3;i++){
		cout<<"K="<<K[i]<<endl;
		for (int j=0;j<4;j++){
			OptionPricing.setK(K[i]);
			OptionPricing.getCosPricing(-bound[j],bound[j]);
		}
		cout<<endl;
	}

	cout<<"BS Formula"<<endl;
	for (int i=0;i<3;i++){
		cout<<"K="<<K[i]<<endl;
			OptionPricing.setK(K[i]);
			OptionPricing.getBSprice();
	}
	return 0;
}