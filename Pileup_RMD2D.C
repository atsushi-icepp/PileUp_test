#include <stdio.h>
#include <iostream>
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#include <TGraph.h>

#define NCONF 8

// Pile up probability is calculated using exp{-deltaT/tau} relation.

#include"gaussian_integrater.h"

void Pileup_RMD2D(){
   // Configurations
   const Double_t width_bundle[NCONF] = {7.,8.,9.,10.,13.,15.,18.,20.};  // Number of fibers in a bundle
   const Double_t deltaT = 50.0e-9;                  // Minimum time difference to distinguish pileup
   const Double_t murate = 9.76E+07;                  // Muon beam rate
   const Double_t mubeam_sigma = 20.39; // sigma of muon beam profile
   const Double_t RMDhit_sigma = 28.13; // sigma of RMD positron hit distribution

   Double_t width_startX=0.; // the strip is located on [width_startX,width_startX+width_bundle], The same is true for Y axis
   Double_t Prob_tot=0.; // total probablity for pileup 
   Double_t Prob_totX = 0.,Prob_totY = 0.; // total probablity for pileup in case the only one strip in the X (Y) direction is implemented
   for (Int_t iconf=0; iconf<NCONF; iconf++) {
      Double_t tauX = 1./murate/gaus_integral(mubeam_sigma,width_startX,width_startX+width_bundle[iconf]);
      Double_t PProbX = 1.-exp(-1.*deltaT/tauX);
      Prob_totX += 2.*PProbX*gaus_integral(RMDhit_sigma,width_startX,width_startX+width_bundle[iconf]);
      Double_t width_startY = 0.;
      for (Int_t jconf=0;jconf<NCONF; jconf++){
         Double_t tauY = 1./murate/gaus_integral(mubeam_sigma,width_startY,width_startY+width_bundle[jconf]);
         Double_t PProbY = 1.-exp(-1.*deltaT/tauY);
         Double_t tauXY = murate*tauX*tauY;
         Double_t PProbXY = 1.-exp(-1.*deltaT/tauXY);
         Double_t PProb = PProbXY + (PProbX-PProbXY)*(PProbY-PProbXY);
         Prob_tot+=4.*PProb*gaus_integral(RMDhit_sigma,width_startX,width_startX+width_bundle[iconf])*gaus_integral(RMDhit_sigma,width_startY,width_startY+width_bundle[jconf]);
         Prob_totY+=4.*PProbY*gaus_integral(RMDhit_sigma,width_startX,width_startX+width_bundle[iconf])*gaus_integral(RMDhit_sigma,width_startY,width_startY+width_bundle[jconf]);
         width_startY += width_bundle[jconf];
      }
      width_startX += width_bundle[iconf];
   } 
   printf("total pile up probability is %lf; simple product is %lf \n",Prob_tot,Prob_totX*Prob_totY);
   printf("total pile up probability in case of only X strip is %lf \n",Prob_totX);
   printf("total pile up probability in case of only Y strip is %lf \n",Prob_totY);
}
