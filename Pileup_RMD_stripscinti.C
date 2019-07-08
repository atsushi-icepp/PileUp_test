#include <stdio.h>
#include <iostream>
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#include <TGraph.h>

#define NCONF 9

// Pile up probability is calculated using exp{-deltaT/tau} relation.

#include"gaussian_integrater.h"

void Pileup_RMD_stripscinti(){
   // Configurations
   const Double_t width_bundle[NCONF] = {4.,4.,4.,4.,6.25,9.,12.25,25.,25.};  // Number of fibers in a bundle
   const Double_t deltaT = 120.0e-9;                  // Minimum time difference to distinguish pileup
   const Double_t murate = 9.76E+07;                  // Muon beam rate
   const Double_t mubeam_sigma = 20.39; // sigma of muon beam profile
   const Double_t RMDhit_sigma = 28.13; // sigma of RMD positron hit distribution

   Double_t width_start=0.; // the fiber is located on [N_start,N_start+Nfib_bundle]
   Double_t Prob_tot=0.;
   for (Int_t iconf=0; iconf<NCONF; iconf++) {
      double tau = 1./murate/gaus_integral(mubeam_sigma,width_start,width_start+width_bundle[iconf]);
      double PProb = 1.-exp(-1.*deltaT/tau);
      Prob_tot+=2.*PProb*gaus_integral(RMDhit_sigma,width_start,width_start+width_bundle[iconf]);
      width_start += width_bundle[iconf];
   } 
   printf("total pile up probability is %lf \n",Prob_tot);
}
