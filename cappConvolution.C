#include <stdio.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TObject.h>
#include <TRandom.h>
#include <TFile.h>
#include <math.h>
#include <TF1Convolution.h>
#include <TF1.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TStopwatch.h>
 
using namespace std;
 
void cappConvolution()
{
   //construction of histogram to fit
   TH1F *h_CappGauss = new TH1F("h_CappGauss","CAPP convoluted by gaussian",100,-1.,6.);
   //TF1 *clapp;
   TF1 *clapp = new TF1("clapp","exp([0]+x*[1])",0.,10.);
   clapp->SetParameter(0,1);
   clapp->SetParameter(1,-0.3);
   /* TCanvas *c1;
   c1 = new TCanvas("c1");
   c1->cd(1);
   clapp->Draw();
   int pippo;
   cin >> pippo;*/
   
   for (int i=0;i<1e6;i++)
   {
      
     Double_t x = clapp->GetRandom(); //gives a alpha of -0.3 in the exp
     x += gRandom->Gaus(0.,0.3);
     h_CappGauss->Fill(x);//probability density function of the addition of two variables is the convolution of 2 dens. functions
   }
 
   TF1Convolution *f_conv = new TF1Convolution("clapp","gaus",-1,6,true);
   f_conv->SetRange(-1.,6.);
   f_conv->SetNofPointsFFT(1000);
   TF1   *f = new TF1("f",*f_conv, 0., 5., f_conv->GetNpar());
   f->SetParameters(1.,-0.3,1,0.,0.3);
   f->FixParameter(2,1);
   //fit
   new TCanvas("c","c",800,1000);

   h_CappGauss->Draw("e");
   h_CappGauss -> Fit("f","","",0,6);
 
}
