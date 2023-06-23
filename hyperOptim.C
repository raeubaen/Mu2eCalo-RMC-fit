#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

using namespace std;

double binw = 0.5; //must be equal to the binw of Resolution histograms from Reso.root
double resoslicew = 5; //in Reso.root è 5
double flmin = -50, flmax=20; // range of legit fluctuactions
int nflbins = (int)((flmax-flmin)/binw);

double ethmin = 30, ethmax = 95; // range of photon energy
int ethbins = (int)((ethmax-ethmin)/binw);
double erecomin = 50, erecomax = 90; // range of reconstructed energy
int erecobins = (int)((erecomax-erecomin)/binw);
double efitmin = 50, efitmax = 90; // fit range
int efitbins = (efitmax - efitmin)/binw;

int npar = 2; // fit parameters (kmax, alpha)

// unfortunately, to use Minuit in thew way we need, global variables are difficult to avoid
TH2D *smearing_pdf;
TF1 *f_eth;
TH1D *eexp_hist;
TFile *resoFile;
double loss_offset;

TH1D *GetErecoHist(const double *par){

  TH1D *ereco_hist = new TH1D("ereco", "ereco", erecobins, erecomin, erecomax);

  f_eth->SetParameters(par[0], par[1]);

  for(int i=0; i<ethbins; i++){
    smearing_pdf->SetAxisRange(ethmin + i*binw, ethmin + (i+1)*binw); //gets slice for current eth bins
    TH1D *tmp = (TH1D*) smearing_pdf->ProjectionY()->Clone(); // projects on ereco direction
    tmp->Scale(f_eth->Eval(ethmin + (i+0.5)*binw)); // multiplies for theoretical pdf
    ereco_hist->Add(tmp); // adds to ereco_hist
    tmp->Delete();
  }

  ereco_hist->Scale(10e6/ereco_hist->Integral()); // normalization
  return ereco_hist;
} // the results must be USED and then FREED, please - I need my RAM



// calculates the histogram of reconstructed energy (from "analytic" smearing) given the parameters
TH1D *GetHist(const double *par){
  //cout << "IN" << endl;
  TH1D *ereco_hist = new TH1D("eexp", "eexp", erecobins, erecomin, erecomax);
  f_eth->SetParameters(par);
  for(int i=0; i<10e6; i++){
    double en = f_eth->GetRandom();
    int eslice = (int)en/resoslicew; //energy slice width for Reso.root
    TH1D *h_res = (TH1D*)resoFile->Get(Form("ALL_cuts/Reso_cutALL_%i", eslice));
    ereco_hist->Fill(en + h_res->GetRandom());
  }

  //ereco_hist->Scale(1/ereco_hist->Integral());
  return ereco_hist;
} // the results must be USED and then FREED, please - I need my RAM


// calculates the chi square from the histograms of MC energy and energy reconstructed by "analytic" smearing, for the current value of parameters
double ChiSquare(const double *par){
   TH1D *ereco_hist = GetErecoHist(par);

   double sum = 0;
   for(int i=0; i<erecobins; i++){
      double x = ereco_hist->GetBinCenter(i);
      //cout << "PARAMETRI: " << par[0] << " " << par[1] << endl;

      //cout << x << endl;

      if(x < efitmin || x > efitmax) continue; // only bins included in fit range

      double y1 = ereco_hist->GetBinContent(i);
      double y2 = eexp_hist->GetBinContent(i);
      double y1err = ereco_hist->GetBinError(i);
      double y2err = eexp_hist->GetBinError(i);

      //cout << y1 << " " << y2 << endl;
/*      double err_sq = (y1err*y1err + y2err*y2err);
      double diff = (y1 - y2) * (y1-y2) / err_sq;
*/

      double log_fact;
      if(y2<30) log_fact = TMath::Log(TMath::Factorial(y2));
      else log_fact = (y2+0.5)*TMath::Log(y2) - y2; // manca una costante

      double diff = y2 * TMath::Log(y1) - y1 - log_fact; // Log(Poisson(obs-MC | th-MC))

      //cout << diff << " ";
      if(!TMath::IsNaN(diff - diff) && y1 !=0) {
        sum -= diff; //+ for chi2, - for nLL
      }
      //cout << sum << endl << endl;
   }
   ereco_hist->Delete(); // to preserve RAM integrity
   double out = sum - loss_offset;
   return out;
}


Double_t clapp(Double_t *xx, Double_t *par)
{
  Float_t x =xx[0];
  double clapp_out;
  if (x < par[0]) clapp_out = (1-2*x/par[0]+2*x/par[0]*x/par[0])*x/par[0]*pow(1-x/par[0], par[1]);
  else clapp_out = 0;
  clapp_out += 0.00045;
  return clapp_out;
}

// main function
double *loss(double ){

  // setting kmax to 90 and alpha to 2
  f_eth->SetParameters(kmax, alpha);

  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Scan");

  // set tolerance , etc...
  minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(100000);  // for GSL
  minimum->SetTolerance(0.1); // very important to tune if the minimization fails - 0.01 is a good value
  minimum->SetPrintLevel(1);

  // create function wrapper for minimizer
  ROOT::Math::Functor f(&ChiSquare, 2);
  minimum->SetFunction(f);

  double step[2] = {kmax*step_perc, alpha*step_perc}; // very important to tune if the minimization fails - 0.005 is a good value
  double variable[2] = { kmax, alpha};    // starting point

  // Set the free variables to be minimized !
  minimum->SetVariable(0,"kmax",variable[0], step[0]);
  minimum->SetVariable(1,"alpha",variable[1], step[1]);

  // do the minimization
  minimum->Minimize();

  double *out = new double[4];

  const double *xs = minimum->X();
  const double *errs = minimum->Errors();

  out[0] = (xs[0] - kmax)/errs[0];
  out[1] = (xs[1] - alpha)/errs[1];
  out[2] = errs[0];
  out[3] = errs[1];

  return out;
}

void hyperOptim(){

  cout.precision(17);

  resoFile = new TFile("Reso.root");

  TFile *sm = new TFile("smearprob.root");
  smearing_pdf = (TH2D*)sm->Get("smearing_pdf");

  f_eth = new TF1("f_eth", clapp, ethmin, ethmax, 2);

  double kmax_min = 90;
  double kmax_max = 91; //91
  double alpha_min = 2;
  double alpha_max = 2.2; //2.2
  int n=1;
  double kmax_w = (kmax_max - kmax_min)/n;
  double alpha_w = (alpha_max - alpha_min)/n;

  double step_perc_min = 0.0005;
  double step_perc_max = 0.05; //0.01
  double mult_step = 1.2589254;
  int nstep = 100;

  double pull_bins[20], step_bins[20], err_bins_k[20] = {}, err_bins_a[20] = {};
  pull_bins[0] = -5;
  step_bins[0] = step_perc_min;

  for(int i=1; i<20; i++){
    pull_bins[i] = pull_bins[i-1] + 0.5;
    step_bins[i] = step_bins[i-1] * mult_step;
    err_bins_k[i] = err_bins_k[i-1] + 0.25;
    err_bins_a[i] = err_bins_a[i-1] + 0.025;
  }

  TH2F *h[4] = {
    new TH2F("pullkmax", "pullkmax", 19, step_bins, 19, pull_bins),
    new TH2F("pullalpha", "pullalpha", 19, step_bins, 19, pull_bins),
    new TH2F("errkmax", "errkmax", 19, step_bins, 19, err_bins_k),
    new TH2F("erralpha", "erralpha", 19, step_bins, 19, err_bins_a)
  };

  for(double kmax=kmax_min; kmax<kmax_max; kmax+=kmax_w){
    for(double alpha=alpha_min; alpha<alpha_max; alpha+=alpha_w){
      eexp_hist = GetHist((double[]){kmax, alpha});
      for(double step_perc=step_perc_min; step_perc<step_perc_max; step_perc*=mult_step){
        double *out = loss(kmax, alpha, step_perc);
        //cout << out[0] << " " << out[1] << " " << out[2] << " " << out[3] << endl;
        for(int i=0; i<4; i++) h[i]->Fill(step_perc, out[i]);
      }
    }
  }

  for(int i=0; i<4; i++){
    TCanvas *c = new TCanvas(Form("c%i", i), Form("c%i", i));
    h[i]->Draw();
  }
}

