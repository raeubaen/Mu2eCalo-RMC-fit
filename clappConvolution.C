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

Double_t resolution=2;
Double_t xfitmin = 10.;
Double_t xfitmax = 100.;
Double_t kmaxrun = 90.;
Double_t funintegral=0;
Double_t correctnorm =1;
Double_t numevfit =0;
int firstloop=0;

Double_t pippofun(Double_t *x, Double_t *par)
{
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t fclapp;
  Double_t xx;
  Double_t sum=0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  int np;
  
  if( firstloop == 0){
    firstloop =1;
    funintegral = 0;
    cout << " FCN  first time call - X: " << x[0] << " par0 " << par[0] << endl;
    cout << " Norm value   " << par[1] << endl;
    cout << " correctnorm: " << correctnorm << endl;
  }
 
  if( par[0] != kmaxrun  ){ // New startup
    //cout << " new FCN round - X " << x[0] << " par-0 " << par[0] << endl;
    //cout << " new FCN round  p1/integ/nev " << par[1] << " " << funintegral << " " << numevfit << endl;
    kmaxrun=par[0];
    
    if(funintegral >0 ){
      if( firstloop==2 ){
	correctnorm=par[1]*numevfit/funintegral;
      }     
      if( firstloop==1  ){     
	correctnorm =  par[1]*numevfit*8/funintegral;
      }
    }     
    // cout << " loop par[1]-correct " << firstloop << " " <<  correctnorm << endl;
    if( firstloop==1) firstloop=2;
    funintegral=0;
  }
  // x is the model value, xlow-xupp the convolution region
  // npoints used for finite convolution integral
  np = 100;
  // Range of convolution integral and step-size
  xlow = x[0] - 5*resolution;
  xupp = x[0] + 5*resolution; 
  step = (xupp-xlow) / np;
 
  // Convolution integral of CLAPP and Gaussian by sum
  float alfa = par[2];
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    double xl = xx/par[0];
    if( xl>0 && xl<1 ){
      fclapp = (1-2*xl+2*xl*xl)*xl*pow((1-xl),alfa);
    }else{
      fclapp=0;
    }
    sum += fclapp * TMath::Gaus(x[0],xx,resolution);
    //cout << " xl " << xl << " " << sum << " fc " << fclapp << " gaus " << TMath::Gaus(x[0],xx,resolution) << endl; 
    
    xx = xupp - (i-.5) * step;
    double xu = xu/par[0];
    if( xu>0 && xu<1){
      fclapp = (1-2*xu+2*xu*xu)*xu*pow((1-xu),alfa);
    }else{
      fclapp=0;
    }
    sum += fclapp * TMath::Gaus(x[0],xx,resolution);
  }
  Double_t  sumscaled;
 
  sumscaled = par[1]*2*step*sum*invsq2pi*0.96/(4.1*resolution);
  funintegral = funintegral+sumscaled;
  //cout << " para x  funintegral : " << par[0] << " " << par[1] << " " << x[0] << " " << funintegral << endl;
  return (sumscaled);
}

void clappConvolution()
{
   //construction of histogram to fit
  TH1F *h_Clapp  = new TH1F("h_Clapp","h_Clapp",100,0.,100.);
  TH1F *h_ClappGauss = new TH1F("h_ClappGauss","CLAPP convoluted by gaussian",110,-10.,100.);
  //TF1 *clapp;
  TF1 *clapp = new TF1("clapp","(1-2*x/[0]+2*x/[0]*x/[0])*x/[0]*(1-x/[0])^2.00",0.,kmaxrun);
  clapp->SetParameter(0,kmaxrun);
  TF1 *clapp1 = new TF1("clapp1","(1-2*x/[0]+2*x/[0]*x/[0])*x/[0]*(1-x/[0])^2.50",0.,90.);
  clapp1->SetParameter(0,kmaxrun);
  
   
   TH1F *resol = new TH1F("resol","resol",1000,-5*resolution,5*resolution);
   TH1F *resolpro = new TH1F("resolpro","resolpro",1000,-5*resolution,5*resolution);
   
   for (int i=0;i<1e5;i++)
   {     
     Double_t x = clapp->GetRandom(); //
     h_Clapp->Fill(x);
     Double_t xres = gRandom->Gaus(0.,resolution);
     x= x+xres;
     resol->Fill(xres);     
     h_ClappGauss->Fill(x);
   }

   //---  Save Histogram to be fit and resolution used --------------------------------------
   TFile *tout;
   tout = new TFile("pippo_out.root","recreate");
   tout->cd();
   h_Clapp->Write();
   h_ClappGauss->Write();
   resol->Write();
   double_t numev = resol->Integral();
   cout << "Numev " << numev << endl;
   int nbin = resol->GetNbinsX();
   Double_t xresmin = resol->GetXaxis()->GetXmin();
   Double_t xresmax = resol->GetXaxis()->GetXmax();
   Double_t stepbin = (xresmax-xresmin)/nbin;
   //cout << " Xmin Xmax Step : " << xresmin << " " << xresmax << " " << stepbin << endl;
   for(int j=0;j<nbin;j++){
    
     double yval = resol->GetBinContent(j)/numev;
     double xval = xresmin + stepbin*j;
     //cout << j << " " << xval  << " " << yval << endl;
     resolpro ->Fill(xval,yval);
     resolpro->SetBinError(j,sqrt(yval)/numev);  
   }
   resolpro->Write();
   tout->Close();
   //===========================================================

   double normaclapp = clapp->Integral(xfitmin,xfitmax);
   cout << " normaclap  " << normaclapp << endl;
   
   int bmin = h_ClappGauss->GetXaxis()->FindBin(xfitmin);
   int bmax = h_ClappGauss->GetXaxis()->FindBin(xfitmax);
   numevfit = h_ClappGauss->Integral(bmin,bmax);
   
   cout << " Numeventtofit: " << numevfit << endl;
   cout << " Fitregion/tot: " << numevfit/numev << endl;
   
   TF1 *ffit = new TF1("ffit",pippofun,xfitmin,xfitmax,3);
   Double_t kmaxstart[3];
   kmaxstart[0]= kmaxrun;
   kmaxstart[1]= numevfit;
   kmaxstart[2]= 2.0;

   ffit->SetParameters(kmaxstart);
   ffit->SetParNames("kmax","Numev","alfa");
   double intfunmio = 0;
   double stepfit =(xfitmax-xfitmin)/100;
   for(int j=0; j<100;j++){
     
     float xcur = xfitmin+stepfit*j;
     float ycur = ffit->Eval(xcur);
     intfunmio = intfunmio + stepfit*ycur;
     // cout << " j " << j << " " << xcur << " " << ycur << " " << intfunmio << endl;

   }
   cout << " myint " << intfunmio << endl;
   double intfun = ffit->Integral(xfitmin,xfitmax);
   
   cout << " IntFun (nocorr): " << intfun << endl;
   double ratio = numevfit/intfun;
   correctnorm =  numevfit*ratio;
   kmaxstart[1]= correctnorm;
   cout << " Ratio (num/fun): " << ratio << " " << correctnorm << endl;
   ffit->SetParameter(1,kmaxstart[1]);
   //ffit->FixParameter(1,kmaxstart[1]*1.105*0.96);
   ffit->FixParameter(2,2.0);
   intfun = ffit->Integral(xfitmin,xfitmax);
   cout << " IntFun (corr) : " << intfun << endl;


   h_ClappGauss->Fit("ffit","RBO");

     
   new TCanvas("c","c",800,1000);

   
   h_ClappGauss->Draw("e");
   ffit->Draw("same");
   //intfun = ffit->Integral(xfitmin,xfitmax);
   //cout << " IntFunFinal " << intfun << endl;
   
   // h_ClappGauss -> Fit("f","","",10,88);
   
 
}

