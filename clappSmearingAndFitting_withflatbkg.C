#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>


using namespace std;

double binw = 0.5; //must be equal to the binw of Resolution histograms from Reso.root

double flmin = -40, flmax=20; // range of legit fluctuactions
int nflbins = (int)((flmax-flmin)/binw); 

double ethmin = 50, ethmax = 95; // range of photon energy
int ethbins = (int)((ethmax-ethmin)/binw); 
double erecomin = 30, erecomax = 110; // range of reconstructed energy 
int erecobins = (int)((erecomax-erecomin)/binw);
double efitmin = 50, efitmax = 110; // fit range
int efitbins = (efitmax - efitmin)/binw;

int npar = 2; // fit parameters (kmax, alpha)

// unfortunately, to use Minuit in thew way we need, global variables are difficult to avoid
TH2F *smearing_pdf_clone;
TF1 *f_eth;
TH1F *eexp_hist;

// calculates the histogram of reconstructed energy (from "analytic" smearing) given the parameters
TH1F *GetErecoHist(const double *par){
  
  TH1F *ereco_hist = new TH1F("ereco", "ereco", erecobins, erecomin, erecomax);

  f_eth->SetParameters(par[0], par[1]);

  for(int i=0; i<ethbins; i++){
    smearing_pdf_clone->SetAxisRange(ethmin + i*binw, ethmin + (i+1)*binw); //gets slice for current eth bins
    TH1F *tmp = (TH1F*) smearing_pdf_clone->ProjectionY()->Clone(); // projects on ereco direction
    double thpdf = f_eth->Eval(ethmin + (i+0.5)*binw);
    tmp->Scale(thpdf); // multiplies for theoretical pdf
    ereco_hist->Add(tmp); // adds to ereco_hist
    tmp->Delete();
  }

  ereco_hist->Scale(1/ereco_hist->Integral()); // normalization
  return ereco_hist;
} // the results must be USED and then FREED, please - I need my RAM


// calculates the chi square from the histograms of MC energy and energy reconstructed by "analytic" smearing, for the current value of parameters
double ChiSquare(const double *par){
   TH1F *ereco_hist = GetErecoHist(par);

   double sum = 0;
   for(int i=0; i<erecobins; i++){

      double x = ereco_hist->GetBinCenter(i);
      if(x < efitmin || x > efitmax) continue; // only bins included in fit range

      double y1 = ereco_hist->GetBinContent(i);
      double y2 = eexp_hist->GetBinContent(i);
      double y1err = ereco_hist->GetBinError(i);
      double y2err = eexp_hist->GetBinError(i);

      double err_sq = (y1err*y1err + y2err*y2err);
      double diff = (y1 - y2) * (y1-y2) / err_sq;
      if(!TMath::IsNaN(diff) && err_sq!=0) sum += diff;
   }
   ereco_hist->Delete(); // to preserve RAM integrity
   return sum/(efitbins - npar - 1);
}

Double_t clapp(Double_t *xx, Double_t *par)
{
   Float_t x =xx[0];
   double clapp_out;
   if (x < par[0]) clapp_out = (1-2*x/par[0]+2*x/par[0]*x/par[0])*x/par[0]*pow(1-x/par[0], par[1]);
   else clapp_out = 0;
   //clapp_out += 0.00045;
   return clapp_out;
}

// main function
void clappSmearingAndFitting_withflatbkg(){

  TFile *resoFile = new TFile("Reso.root");

  int retrieveSmearingProb = 0;
  TH2F *smearing_pdf;

  if(retrieveSmearingProb == 0){
    cout << "Evaluating smearing matrix" << endl;
    // smearing matrix
    smearing_pdf = new TH2F("smearing_pdf", "Smearing PDF", ethbins, ethmin, ethmax, erecobins, erecomin, erecomax);
    for(int i=0; i<ethbins; i++){
      double eth = ethmin + i*binw;
      int eslice = (int)eth/5; //energy slice width for Reso.root

      TH1F *h_res = (TH1F*)resoFile->Get(Form("ALL_cuts/Reso_cutALL_%i", eslice));
      h_res->SetAxisRange(flmin, flmax); // adjusting the range to speed up
      h_res->Scale(1/h_res->Integral()); // normalization
      for(int j=0; j<nflbins; j++){

        double fluct = flmin + j*binw;
        double ereco = eth + fluct;
        smearing_pdf->SetBinContent(i, (int)((ereco - erecomin)/binw), h_res->GetBinContent(h_res->FindBin(fluct)));
        smearing_pdf->SetBinError(i, (int)((ereco - erecomin)/binw), h_res->GetBinError(h_res->FindBin(fluct)));
      }
    }

    TFile *sm = new TFile("smearprob.root", "recreate");
    sm->cd();
    smearing_pdf->Write();
  }

  else{
    TFile *sm = new TFile("smearprob.root");
    smearing_pdf = (TH2F*)sm->Get("smearing_pdf");
  }

  TCanvas *c1 = new TCanvas("c1", "Smearing Matrix");
  c1->cd();
  smearing_pdf->Draw("zcol");

  smearing_pdf_clone = (TH2F*)smearing_pdf->Clone();

  double kmax = 90, alpha=2;

  TCanvas *c2 = new TCanvas("c2", "closure approximation function");
  // clapp function, having kmax and alpha as parameters
  f_eth = new TF1("f_eth", clapp, ethmin, ethmax, 2);
  // setting kmax to 90 and alpha to 2
  f_eth->SetParameters(kmax, alpha);
  TF1 *f_eth_clone = (TF1*) f_eth->Clone();
  f_eth_clone->Draw();

  cout << "Simulating detector response (Montecarlo)" << endl;
  // MC extraction
  eexp_hist = new TH1F("eexp", "eexp", erecobins, erecomin, erecomax);
  for(int i=0; i<1e6; i++){
    double en = f_eth->GetRandom();
    int eslice = (int)en/5; //energy slice width for Reso.root
    if(eslice < 0 || eslice>19) continue;
    TH1F *h_res = (TH1F*)resoFile->Get(Form("ALL_cuts/Reso_cutALL_%i", eslice));
    eexp_hist->Fill(en + h_res->GetRandom());
  }
  eexp_hist->Scale(1/eexp_hist->Integral()); // normalization


  gErrorIgnoreLevel = kError;
  // construction minuit minimizer
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

  // set tolerance , etc...
  minimum->SetMaxFunctionCalls(10000000); // for Minuit/Minuit2
  minimum->SetMaxIterations(100000);  // for GSL
  minimum->SetTolerance(0.01); // very important to tune if the minimization fails - 0.01 is a good value
  minimum->SetPrintLevel(1);

  // create function wrapper for minimizer
  ROOT::Math::Functor f(&ChiSquare, 2);
  minimum->SetFunction(f);

  double step[2] = {0.005,0.005}; // very important to tune if the minimization fails - 0.005 is a good value
  double variable[2] = { kmax, alpha};    // starting point

  // Set the free variables to be minimized !
  minimum->SetVariable(0,"kmax",variable[0], step[0]);
  minimum->SetVariable(1,"alpha",variable[1], step[1]);

  cout << endl << "Fitting in range " << efitmin << ", " << efitmax << endl << endl;

  // do the minimization
  minimum->Minimize();

  gErrorIgnoreLevel = kWarning;

  const double *xs = minimum->X();
  std::cout << "Minimum: Chi-square(" << xs[0] << "," << xs[1] << "): "
            << minimum->MinValue()  << std::endl;

  c2->cd();
  f_eth->SetParameters(xs);
  f_eth->SetLineColor(kBlack);
  f_eth->Draw("same");

  TH1F *ereco_hist = GetErecoHist(xs);

  TCanvas *c3 = new TCanvas("c3", "Energy histograms (also outside the fit range)");

  ereco_hist->SetLineColor(kRed);
  eexp_hist->Draw();
  ereco_hist->Draw("same");

  TLegend *l = new TLegend();
  l->AddEntry(ereco_hist, "convolution");
  l->AddEntry(eexp_hist, "MC/data");
  l->Draw();

  TCanvas *c4 = new TCanvas("c4", "Energy histograms (in the fit range)");

  TH1F *ereco_hist_clone = (TH1F*) ereco_hist->Clone();
  ereco_hist_clone->SetAxisRange(efitmin, efitmax);
  ereco_hist_clone->SetLineColor(kRed);
  eexp_hist->Draw();
  ereco_hist_clone->Draw("same");
  l->Draw();

  TString outFileName = "out.root";
  TFile *outFile = new TFile(outFileName, "recreate");
  outFile->cd();
  cout << endl << "Output in " << outFileName.Data() << endl << endl;
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  outFile->Close();
}
