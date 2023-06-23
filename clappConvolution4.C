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
int  numevgen= 2000000;
int flag=1; // if 1 Reso.root is used, otherwise a gaussian 
Double_t resolution=5.0; 
Double_t xfitmin = 40.;
Double_t xfitmax = 100.;
Double_t kmaxrun = 90.;
Double_t funintegral=0;
Double_t correctnorm =1;
Double_t numevfit =0;
int firstloop=0;
Double_t ProbReso[20][300];
Double_t XminReso[20];
TH1F* hh[20];  

int MakeSmearingProb(){ //converts the histograms to vectors containing probabilities
  int status = 0;
  cout << " reading Resolution file Reso.root" << endl;
  TFile *f1 = new TFile("Reso.root");
  if( f1!= NULL ){ // if File is found status= success
    status = 1;
    f1->ls();
    int jgam;
    for (int jrun=0; jrun<20; jrun++){
      if( jrun<2) {
	      jgam=2;
      } else{
	      jgam=jrun;
      }
      hh[jrun] = (TH1F*)f1->Get(Form("ALL_cuts/Reso_cutALL_%d",jgam));
      cout << " jrun " << jrun << " jgam .. " << jgam << " .. loaded " << endl;      
      //pgam = jgam*5;

      double NevTot = hh[jrun]->Integral();
      int NumBin = hh[jrun]->GetNbinsX();
      float xmin = hh[jrun]->GetXaxis()->GetXmin();
      float xmax = hh[jrun]->GetXaxis()->GetXmax();
      XminReso[jrun] = xmin;
      float step = (xmax-xmin)/NumBin;     
      double xval[300];
      double yval[300];
      //--------------------------------------------- 
      for(int j=0;j<NumBin;j++){
        if( NevTot ==0){
          ProbReso[jrun][j] = 0;
        } else{
          xval[j] = hh[jrun]->GetBinCenter(j+1);
          yval[j] = hh[jrun]->GetBinContent(j+1);
          ProbReso[jrun][j]= yval[j]/NevTot;
          cout << " jrun j " << jrun << " " << j << " Prob " << ProbReso[jrun][j] << endl;
        }
      } // Loop over 300 reconstructed energy bins
    }   // Loop over all 20 bands of true energies
  }
  return status;
}

//============================================================
Double_t fun(Double_t *x, Double_t *par)
{
  //int flag=0; // 0=gaussian, 1 -reso 
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

  if( flag==1){
    np=120;
    step = 0.5;
    xlow = x[0]-50;
    xupp = x[0]+10;
    float alfa = par[2];
    //-------------------------------------------
    // xx is the model, x[0] is the reconstructed
    //-------------------------------------------
    double distance;
    int  igam, ibinreso;
    
    for(int i=1; i <120; i++){
      xx = xlow + (i-1)*step + step/2;
      igam = xx/5;
      float distance = x[0]-xx;
      float xmin = XminReso[igam];
      if( distance-xmin > 0){
        ibinreso = (distance-xmin)/step;
        double smearprob = ProbReso[igam][ibinreso];
          
        double xl = xx/par[0];
        //cout << " x[0] xx " << x[0] << " " << xx << " "  << " dis ib " << distance << " " << ibinreso << "  prob " << smearprob << endl;
        if( xl>0 && xl<1 ){
          fclapp = (1-2*xl+2*xl*xl)*xl*pow((1-xl),alfa);	
        } else{
          fclapp=0;
        }
        sum += fclapp*smearprob;
      }
    }     
    //cout << " sum " << sum << " fclap " << fclapp << " xrec " << x[0] << endl;
    
  }else if( flag==0){ 
    np = 100;
    //---------------------------------------------
    // Range of convolution integral and step-size
    //---------------------------------------------
    xlow = x[0] - 5*resolution;
    xupp = x[0] + 5*resolution; 
    step = (xupp-xlow) / np;
    
    float alfa = par[2];
    for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      double xl = xx/par[0];
      if( xl>0 && xl<1 ){
	      fclapp = (1-2*xl+2*xl*xl)*xl*pow((1-xl),alfa);
      } else{
	      fclapp=0; 
      }
      sum += fclapp * TMath::Gaus(x[0],xx,resolution);
      
      xx = xupp - (i-.5) * step;
      double xu = xx/par[0];  // bug found .. xu = xl/par[0] now xu= xx/par[0]
      if( xu>0 && xu<1){
	      fclapp = (1-2*xu+2*xu*xu)*xu*pow((1-xu),alfa);
      } else{
	      fclapp=0;
      }      
      sum += fclapp * TMath::Gaus(x[0],xx,resolution);
    }
  }
  //-----------------------------------------------
  
  Double_t  sumscaled;
  if( flag==1){
    sumscaled = par[1]*step*sum;
  }else if(flag==0){
    //0.2222 is the clapp integral for 0.1 resolution ..
    sumscaled = par[1]*step*sum*invsq2pi/(4.55272*resolution); // 4.55272 è 1/0.222? non esattamente -- ???...
  }
  funintegral = funintegral+sumscaled;
  //cout << " para x  funintegral : " << par[0] << " " << par[1] << " " << x[0] << " " << funintegral << endl;
  return (sumscaled);
}

void clappConvolution4(){

  int probstatus;
  
  probstatus = MakeSmearingProb();
  if( probstatus == 0) cout << "nothing loaded" << endl;
  if( probstatus == 1) cout << " Resolution loaded " << endl;
  
  //
  //construction of histogram to fit
  // ===================================
  TH1F *h_Clapp  = new TH1F("h_Clapp","h_Clapp",100,0.,100.);
  TH1F *h_ClappGauss = new TH1F("h_ClappGauss","CLAPP convoluted by gaussian",110,-10.,100.);
  TH1F *h_ClappReso  = new TH1F("h_ClappReso","CLAPP convoluted by resolution",110,-10.,100.);
  
  TF1 *clapp = new TF1("clapp","(1-2*x/[0]+2*x/[0]*x/[0])*x/[0]*(1-x/[0])^2.00",0.,kmaxrun);
  clapp->SetParameter(0, kmaxrun);
   
   TH1F *resol = new TH1F("resol","resol",1000,-5*resolution,5*resolution);
   TH1F *resolpro = new TH1F("resolpro","resolpro",1000,-5*resolution,5*resolution);
   
   for (int i=0;i<numevgen;i++) // 100 keventi generati
   {     
     Double_t x = clapp->GetRandom(); //
     h_Clapp->Fill(x);
     Double_t pgam = x;   
     Double_t xres = gRandom->Gaus(0.,resolution);
     x= x+xres;
     resol->Fill(xres);     
     h_ClappGauss->Fill(x);
     int jgam = pgam/5;
     Double_t xrescalo = hh[jgam]->GetRandom();
     Double_t erec = pgam+xrescalo;
     h_ClappReso->Fill(erec);
     
   }

   //---  Save Histogram to be fit and resolution used --------------------------------------
   TFile *tout;
   tout = new TFile("out.root","recreate");
   tout->cd();
   h_Clapp->Write();
   h_ClappGauss->Write();
   h_ClappReso->Write();
   resol->Write();
   double_t numev = resol->Integral();
   cout << "Numev generated " << numev << endl;
   
   int nbin = resol->GetNbinsX();
   Double_t xresmin = resol->GetXaxis()->GetXmin();
   Double_t xresmax = resol->GetXaxis()->GetXmax();
   Double_t stepbin = (xresmax-xresmin)/nbin;

   for(int j=0;j<nbin;j++){
     double yval = resol->GetBinContent(j)/numev;
     double xval = xresmin + stepbin*j;
     //cout << j << " " << xval  << " " << yval << endl;
     resolpro ->Fill(xval,yval);
     resolpro->SetBinError(j,sqrt(yval)/numev);  
   }
   resolpro->Write();
   tout->Close();
   //
   //return;
   //===========================================================

   if( flag == 0){
     double normaclapp = clapp->Integral(xfitmin,xfitmax);
     cout << " Normalization clap-no-smeared  " << normaclapp << endl;
     
     int bmin = h_ClappGauss->GetXaxis()->FindBin(xfitmin);
     int bmax = h_ClappGauss->GetXaxis()->FindBin(xfitmax);
     numevfit = h_ClappGauss->Integral(bmin,bmax);
     
     cout << " Numeventtofit (smeared): " << numevfit << endl;
     cout << " Fitregion/tot (smeared): " << numevfit/numev << endl;
     
     TF1 *ffit = new TF1("ffit",fun,xfitmin,xfitmax,3);
     Double_t kmaxstart[3];
     kmaxstart[0]= kmaxrun;
     kmaxstart[1]= numevfit;
     kmaxstart[2]= 2.0;
     
     ffit->SetParameters(kmaxstart);
     ffit->SetParNames("kmax","Numev","alfa");
     double intfunmio = 0;
     int nstepmio = 1000;
     double stepfit =(xfitmax-xfitmin)/nstepmio;
     for(int j=0; j<nstepmio;j++){     
       float xcur = xfitmin+stepfit*j;
       float ycur = ffit->Eval(xcur);
       intfunmio = intfunmio + stepfit*ycur;
       // cout << " j " << j << " " << xcur << " " << ycur << " " << intfunmio << endl;
     }
     
     cout << " Manual Integration of smeared function before fit " << intfunmio << endl;
     double ratio = numevfit/intfunmio;
     correctnorm =  numevfit*ratio;
     
     kmaxstart[1]= correctnorm;
     cout << " Ratio : " << ratio << " " << correctnorm << endl;
     ffit->SetParameter(1,kmaxstart[1]);
     
     ffit->FixParameter(1,kmaxstart[1]);
     h_ClappGauss->Fit("ffit","RBO");
     
     new TCanvas("c","c",800,1000);
     
   
     h_ClappGauss->Draw("e");
     ffit->Draw("same");
     h_ClappGauss->Write();
     tout->Close();
     //intfun = ffit->Integral(xfitmin,xfitmax);
     //cout << " IntFunFinal " << intfun << endl;
   
     // h_ClappGauss -> Fit("f","","",10,88);
     
     intfunmio = 0;
     double intfunhigh =0;
     stepfit =(xfitmax-xfitmin)/nstepmio;
     for(int j=0; j<nstepmio;j++){     
       float xcur = xfitmin+stepfit*j;
       float ycur = ffit->Eval(xcur);
       intfunmio = intfunmio + stepfit*ycur;
       if(xcur>=57) intfunhigh = intfunhigh+stepfit*ycur;
       // cout << " j " << j << " " << xcur << " " << ycur << " " << intfunmio << endl;
     }
     cout << " Manual Integration of smeared function after fit " << intfunmio << endl;
     cout << " Manual Integration between 57 up " << intfunhigh << endl;
     
   }else if(flag==1){

     double normaclapp = clapp->Integral(xfitmin,xfitmax);
     cout << " Normalization clap-no-smeared  " << normaclapp << endl;
     
     int bmin = h_ClappReso->GetXaxis()->FindBin(xfitmin);
     int bmax = h_ClappReso->GetXaxis()->FindBin(xfitmax);
     numevfit = h_ClappReso->Integral(bmin,bmax);
     
     cout << " Numeventtofit (smeared): " << numevfit << endl;
     cout << " Fitregion/tot (smeared): " << numevfit/numev << endl;
     
     TF1 *ffit = new TF1("ffit", fun,xfitmin,xfitmax,3);
     Double_t kmaxstart[3];
     kmaxstart[0]= kmaxrun;
     kmaxstart[1]= numevfit;
     kmaxstart[2]= 2.0; // alfa
    
     
     ffit->SetParameters(kmaxstart);
     ffit->SetParNames("kmax","Numev","alfa");
     double intfunmio = 0;
     int nstepmio = 1000;
     double stepfit =(xfitmax-xfitmin)/nstepmio;
     for(int j=0; j<nstepmio;j++){     
       float xcur = xfitmin+stepfit*j;
       float ycur = ffit->Eval(xcur);
       intfunmio = intfunmio + stepfit*ycur;
       //cout << " j " << j << " " << xcur << " " << ycur << " " << intfunmio << endl;
     }
    
     cout << " Manual Integration of smeared function before fit " << intfunmio << endl;
     double intfun = ffit->Integral(xfitmin,xfitmax);
     cout << " Root  Integration of function all range: " << intfun << endl;
     double ratio = numevfit/intfunmio;
     correctnorm =  numevfit*ratio;
     
     kmaxstart[1]= correctnorm;
     cout << " Ratio : " << ratio << " " << correctnorm << endl;
     ffit->SetParameter(1,kmaxstart[1]);   
     ffit->FixParameter(1,kmaxstart[1]);
     
     //ffit->FixParameter(2,2.0);
     //intfun = ffit->Integral(xfitmin,xfitmax);
     //cout << " IntFun (corr) : " << intfun << endl;
     h_ClappReso->Fit("ffit","RBO");
     
     new TCanvas("c","c",800,1000);
     
   
     h_ClappReso->Draw("e");
     ffit->Draw("same");
     h_ClappReso->Write();
     tout->Close();
     //intfun = ffit->Integral(xfitmin,xfitmax);
     //cout << " IntFunFinal " << intfun << endl;
     
     // h_ClappGauss -> Fit("f","","",10,88);
     
     intfunmio = 0;
     double intfunhigh =0;
     stepfit =(xfitmax-xfitmin)/nstepmio;
     for(int j=0; j<nstepmio;j++){     
       float xcur = xfitmin+stepfit*j;
       float ycur = ffit->Eval(xcur);
       intfunmio = intfunmio + stepfit*ycur;
       if(xcur>=57) intfunhigh = intfunhigh+stepfit*ycur;
       // cout << " j " << j << " " << xcur << " " << ycur << " " << intfunmio << endl;
     }
     cout << " Manual Integration of smeared function after fit " << intfunmio << endl;
     cout << " Manual Integration between 57 up " << intfunhigh << endl;
   }
}

