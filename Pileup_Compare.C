#include <stdio.h>
#include <iostream>
#include <TH1.h>
#include <TMath.h>
#include <TFile.h>
#include <TRandom.h>
#include <TGraph.h>

#define NCONF 7

// Pile up probability is calculated using MC and exp{-deltaT/tau} relation.
// They are plotted to the same graph and can be compared.

#include"gaussian_integrater.h"

void Pileup_Compare(Int_t nevents=10000)
{

   // Configurations
   Int_t Nfib_bundle[NCONF] = {4,9,16,25,36,49,64};  // Number of fibers in a bundle
   Int_t Nfib_all = 392;                        // Number of all fibers
   Double_t deltaT = 120.0e-9;                  // Minimum time difference to distinguish pileup
   Double_t murate = 9.76E+07;                  // Muon beam rate

   // Number of pileup events at each bundle
   Double_t Npileup[NCONF][Nfib_all];
   for (Int_t iconf=0; iconf<NCONF; iconf++) {
      for(Int_t iedge=0; iedge<Nfib_all; iedge++){
         Npileup[iconf][iedge] = 0.0;
      }
   }

   // Number of RMD positrons at each bundle (for no
   Double_t N_RMD[NCONF][Nfib_all];
    for (Int_t iconf=0; iconf<NCONF; iconf++) {
      for(Int_t iedge=0; iedge<Nfib_all; iedge++){
         N_RMD[iconf][iedge] = 0.0;
      }
   }   

   // Loop for RMD events
   for (Int_t ievent = 0; ievent < nevents; ievent++) {

      if (ievent % 1000 == 0) {
         printf("Event: %d\r", ievent);
         fflush(stdout);
      }

      // RMD hit timing and position
      Double_t tRMD = 0;
      Double_t zRMD = fabs(gRandom->Gaus(0, 116.5));

      // Muon beam hit timing and position
      Int_t Nmu = gRandom->Poisson(500e-9 * murate/2.); // Only upper half is taken into account
      std::vector<Double_t> tMu;
      std::vector<Double_t> zMu;
      for (Int_t imu = 0; imu < Nmu; imu++) {
         tMu.push_back(gRandom->Rndm() * 500e-9 - 250e-9); // center of tMu is set to zero
         zMu.push_back(fabs(gRandom->Gaus(0, 85.0)));
      }
                  
      // Loop for bundle configurations
      for (Int_t iconf = 0; iconf < NCONF; iconf++) {

         // Identify muon pileup events
         for(Int_t iedge=0; iedge < Nfib_all; ++iedge){
            if(iedge <= zRMD && zRMD <= iedge+Nfib_bundle[iconf] && zRMD < Nfib_all){
               
               N_RMD[iconf][iedge] += 1.;
               for (Int_t imu = 0; imu < Nmu; imu++) {
                  if(iedge <= zMu[imu] && zMu[imu] < iedge+Nfib_bundle[iconf] && TMath::Abs(tRMD-tMu[imu])<deltaT/2.0){
                     Npileup[iconf][iedge] += 1;
                     break;
                  }
               }
            }
         }
         
      }
    
   }
   printf("\n");
   TCanvas *canvEff;
   TGraph *hoge[NCONF];
   TLegend *legend = new TLegend( 0.6, 0.48, 0.99, 0.93);
   canvEff = new TCanvas("RMDinefficiency", "RMD pileup ineffciency",850, 650);
   for (Int_t iconf=0; iconf<NCONF; iconf++) {
      hoge[iconf] = new TGraph();
      hoge[iconf]->SetTitle("Pileup inefficiency vs distance from the origin; distance [cm];Pileup probability");
      hoge[iconf]->SetMaximum(1.);hoge[iconf]->SetMinimum(0.);
      hoge[iconf]->SetMarkerColor(iconf+1);
      hoge[iconf]->SetMarkerSize(0.5);
      hoge[iconf]->SetMarkerStyle(22);
      for(Int_t iedge=0; iedge<Nfib_all; iedge++){
         hoge[iconf]->SetPoint(iedge,iedge*0.025,Npileup[iconf][iedge]/N_RMD[iconf][iedge]);
      }
      if (iconf==0) hoge[iconf]->Draw("ap");
      if (iconf!=0) hoge[iconf]->Draw("p same");
      legend->AddEntry( hoge[iconf], Form("0.25mm*%d strip width(MC)",Nfib_bundle[iconf]) , "p") ;
   } 

// ==== calculation using the equation of exp{-deltaT/tau} ========
   TGraph *foo[NCONF];
   for (Int_t iconf=0; iconf<NCONF; iconf++) {
      foo[iconf] = new TGraph();
      foo[iconf]->SetMaximum(1.);foo[iconf]->SetMinimum(0.);
      foo[iconf]->SetLineColor(iconf+NCONF+1);
      foo[iconf]->SetMarkerSize(0.5);
      foo[iconf]->SetMarkerStyle(8);
      for(Int_t iedge=0; iedge<Nfib_all; iedge++){
         double tau = 1./murate/gaus_integral(85.0,(double)iedge,(double)iedge+Nfib_bundle[iconf]);
         double PEff = 1.-exp(-1.*deltaT/tau);
         foo[iconf]->SetPoint(iedge,iedge*0.025,PEff);
      }
      foo[iconf]->Draw("same");
      legend->AddEntry(foo[iconf], Form("0.25mm*%d strip width(Analytical)",Nfib_bundle[iconf]) , "l") ; 
   } 
   legend->Draw("same");
   canvEff->SaveAs("inefficiency_plot.pdf");
}
