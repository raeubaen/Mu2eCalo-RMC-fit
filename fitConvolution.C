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
 
void fitConvolution(float xmin, float xmax, float xfit1, float xfit2)
{
   //construction of histogram to fit
  //TH1F *h_Exp =new TH1F("h_Exp","h_Exp",100,0.,5.);
  TH1F *h_ExpGauss = new TH1F("h_ExpGauss","Exponential convoluted by gaussian",100,xmin,xmax);
  for (int i=0;i<1e6;i++)
    {
      Double_t x = gRandom->Exp(+1/0.3);//gives a alpha of -0.3 in the exp
      //h_Exp->Fill(x);
      x += gRandom->Gaus(0.,0.3);
      h_ExpGauss->Fill(x);//probability density function of the addition of two variables is the convolution of 2 dens. functions
   }
 
   TF1Convolution *f_conv = new TF1Convolution("expo","gaus",xfit1,xfit2,true);
   f_conv->SetRange(-1,6.);
   f_conv->SetNofPointsFFT(10000);
   TF1   *f = new TF1("f",*f_conv, 0., 5., f_conv->GetNpar());
   f->SetParameters(+1.,-0.3,0.,0.3);
 
   //fit
   new TCanvas("c","c",800,1000);
   // h_Exp->Draw();
   
   h_ExpGauss->Draw("e");
   h_ExpGauss -> Fit("f"," ","",xfit1,xfit2);
   //h_Exp->Draw("same"); 
}
