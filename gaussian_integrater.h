#include <stdio.h>
#include <iostream>
#include <TMath.h>
#include <TRandom.h>

// This header file is used to calculate the integral of normalized gauss distribution with mu=0, while sigma can be input by users.

const int n_Simpson = 7; // This determines the precision of the integral of the Simpson method.
const double SQRT1_2PI = 1./TMath::Sqrt(2.*TMath::Pi());
double normalized_gaus(double x){return SQRT1_2PI*exp(-x*x/2);}

double gaus_integral(const double sigma,const double a,const double b){
// Integral is calculated with the functional values at 2^{n_Simpson} points; i.e.  f(a)+4f(a+h)+f(a+2h) is calculated at 2^{n_Simpson-1} districts, whete  h = (b-a)/2^{n_Simpson -1}
   int i;
   int N = 1;
   for (i=0;i<n_Simpson;i++) N*=2;
   double h = (b-a)/N;
   double integral = 0.;
   for (i=1;i<N;i+=2) integral += 4.*normalized_gaus((a+h*i)/sigma);
   for (i=2;i<N;i+=2) integral += 2.*normalized_gaus((a+h*i)/sigma);
   integral += normalized_gaus(a/sigma)+normalized_gaus(b/sigma);
   integral *= h/(3.*sigma);
   return integral;
}
