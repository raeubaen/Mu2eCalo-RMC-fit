#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include <iostream>

using namespace std;

double binw = 0.5; //must be equal to the binw of Resolution histograms from Reso.root
double resoslicew = 0.5; //in Reso.root è 5
double flmin = -60, flmax=20; // range of legit fluctuactions
int nflbins = (int)((flmax-flmin)/binw); 

double ethmin = 40, ethmax = 150; // range of photon energy
int ethbins = (int)((ethmax-ethmin)/binw); 
double erecomin = 0, erecomax = 170; // range of reconstructed energy 
int erecobins = (int)((erecomax-erecomin)/binw);
double efitmin = 60, efitmax = 110; // fit range
int efitbins = (efitmax - efitmin)/binw;

int npar = 2; // fit parameters (kmax, alpha)

// unfortunately, to use Minuit in thew way we need, global variables are difficult to avoid
TH2D *smearing_pdf_clone;
TF1 *f_eth;
TH1D *eexp_hist;
TFile *resoFile;
TH2D *chi2;
TGraph *path;
int k;
double loss_offset;

TH1D *GetErecoHist(const double *par){

  TH1D *ereco_hist = new TH1D("ereco", "ereco", erecobins, erecomin, erecomax);

  f_eth->SetParameters(par[0], par[1]);

  for(int i=0; i<ethbins; i++){
    smearing_pdf_clone->SetAxisRange(ethmin + i*binw, ethmin + (i+1)*binw); //gets slice for current eth bins
    TH1D *tmp = (TH1D*) smearing_pdf_clone->ProjectionY()->Clone(); // projects on ereco direction
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

   //cout << k << endl;
   //k++;
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
   chi2->SetBinContent(chi2->FindBin(par[0], par[1]), out);
   path->AddPoint(par[0], par[1]);
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
void flatbkg_notbad(){

  cout.precision(17);

  k=0;
  resoFile = new TFile("../outfitReso.root");

  int retrieveSmearingProb = 0;
  TH2D *smearing_pdf;

  if(retrieveSmearingProb == 0){
    cout << "Evaluating smearing matrix" << endl;
    // smearing matrix
    smearing_pdf = new TH2D("smearing_pdf", "Smearing PDF", ethbins, ethmin, ethmax, erecobins, erecomin, erecomax);
    for(int i=0; i<ethbins; i++){
      double eth = ethmin + i*binw;
      int eslice = (int)(eth/resoslicew); //energy slice width for Reso.root
      TH1D *h_res = (TH1D*)resoFile->Get(Form("ALL_cuts/Reso_cutALL_%i", eslice));
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
    smearing_pdf = (TH2D*)sm->Get("smearing_pdf");
  }

  TCanvas *c1 = new TCanvas("c1", "Smearing Matrix");
  c1->cd();
  smearing_pdf->Draw("zcol");

  smearing_pdf_clone = (TH2D*)smearing_pdf->Clone();

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

  double par[2] = {kmax, alpha};
  eexp_hist = GetHist(par);

  double lim[2][2] = {{89.8, 90.2}, {1.98, 2.02}};
  double nbins = 40;
  double w[2];
  for(int i=0; i<2; i++) w[i] = (lim[i][1]-lim[i][0])/nbins;

  chi2 = new TH2D("h", "h", nbins, lim[0][0], lim[0][1], nbins, lim[1][0], lim[1][1]); //era bello 88.9-90-1 1.99-2.01
  path = new TGraph();

  loss_offset = 0;
  cout << "Bruteforcing" << endl;

  double brutal_min = ChiSquare((double[]){(lim[0][0]+lim[0][1])/2, (lim[1][0]+lim[1][1])/2});
  double brutal_min_par[2];
  for(double km=lim[0][0]; km<lim[0][1]; km+=w[0]){
    for(double a=lim[1][0]; a<lim[1][1]; a+=w[1]){
      double loss = ChiSquare((double[]){km, a});
      if(0 < loss && loss < brutal_min){ brutal_min = loss; brutal_min_par[0]=km; brutal_min_par[1] = a; };
    }
  }

  loss_offset = brutal_min;

  cout << "Brute-force minimum: " << brutal_min << endl;
  cout << "Minimum parameter: " << brutal_min_par[0] << " " << brutal_min_par[1] << endl;

  TH2D *chi2_clone = (TH2D*)chi2->Clone();

  // construction minuit minimizer
  ROOT::Math::Minimizer* minimum = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Hesse");

  // set tolerance , etc...
  minimum->SetMaxFunctionCalls(100); // for Minuit/Minuit2
  minimum->SetMaxIterations(100);  // for GSL
  minimum->SetTolerance(0.1); // very important to tune if the minimization fails - 0.01 is a good value
  minimum->SetPrintLevel(1);
  minimum->SetPrecision(1e-14);

  // create function wrapper for minimizer
  ROOT::Math::Functor f(&ChiSquare, 2);
  minimum->SetFunction(f);

  double step[2] = {0.1, 0.01}; // very important to tune if the minimization fails - 0.005 is a good value
  double variable[2] = { kmax+0.05, alpha+0.01};    // starting point

  // Set the free variables to be minimized !
  minimum->SetVariable(0,"kmax",variable[0], step[0]);
  minimum->SetVariable(1,"alpha",variable[1], step[1]);
  minimum->SetVariableLimits(0, 88, 92);
  minimum->SetVariableLimits(1, 1.8, 2.2);

  cout << endl << "Fitting in range " << efitmin << ", " << efitmax << endl << endl;

  path->Delete();
  path = new TGraph();

  // do the minimization
  minimum->Minimize();

  const double *xs = minimum->X();

  double fcn = ChiSquare(xs);
  cout <<endl << endl << fcn << endl << endl;

  loss_offset = - fcn + brutal_min;

  TGraph *path_clone = (TGraph*)path->Clone();

  double x[100], y[100];
  unsigned int nstep = 100;
  minimum->Scan(0, nstep, x, y);
  TGraph *g0 = new TGraphErrors(100, x, y);
  TCanvas *c6 = new TCanvas("c6", "c");
  g0->Draw();
  minimum->Scan(1, nstep, x, y);
  TGraph *g1 = new TGraphErrors(100, x, y);
  TCanvas *c7 = new TCanvas("c7", "c");
  g1->Draw();

  double xx[100], yy[100];
  TCanvas *c8 = new TCanvas("c8", "contour");
  minimum->Contour(0, 1, nstep, xx, yy);
  TGraph *g2 = new TGraph(100, xx, yy);
  g2->Draw();

  gErrorIgnoreLevel = kWarning;

  const double *errs = minimum->Errors();
  std::cout << "Minimum: Chi-square(" << xs[0] << "," << xs[1] << "): "
            << minimum->MinValue()  << std::endl;

  cout << "kmax pull: " << (xs[0] - kmax)/errs[0] << endl;
  cout << "alpha pull: " << (xs[1] - alpha)/errs[1] << endl;

  cout << "Correlation Coefficient between kmax and alpha: " << minimum->Correlation(0, 1) << endl;

  cout << "FCN at kmax+kmax_err, alpha: " << ChiSquare((double[]){xs[0]+errs[0], xs[1]}) << endl;;
  cout << "FCN at kmax, alpha+alpha_err: " << ChiSquare((double[]){xs[0], xs[1]+errs[1]}) << endl;

  c2->cd();
  f_eth->SetParameters(xs);
  f_eth->SetLineColor(kBlack);
  f_eth->Draw("same");

  TH1D *ereco_hist = GetHist(xs);

  TCanvas *c3 = new TCanvas("c3", "Energy histograms (also outside the fit range)");

  ereco_hist->SetLineColor(kRed);
  eexp_hist->Draw();
  ereco_hist->Draw("same");

  TLegend *l = new TLegend();
  l->AddEntry(ereco_hist, "convolution");
  l->AddEntry(eexp_hist, "MC/data");
  l->Draw();

  TCanvas *c4 = new TCanvas("c4", "Energy histograms (in the fit range)");

  TH1D *ereco_hist_clone = (TH1D*) ereco_hist->Clone();
  ereco_hist_clone->SetAxisRange(efitmin, efitmax);
  ereco_hist_clone->SetLineColor(kRed);
  eexp_hist->Draw();
  ereco_hist_clone->Draw("same");
  l->Draw();

  double chi2min = chi2_clone->GetMinimum();

  for(int i=0; i<nbins; i++){
    for(int j=0; j<nbins; j++){
      chi2_clone->SetBinContent(i, j, chi2_clone->GetBinContent(i, j) - chi2min + 0.1);
    }
  }

  TCanvas *c5 = new TCanvas("c5", "");
  c5->cd();
  chi2_clone->Draw("zcol");
  path_clone->SetMarkerStyle(22);
  path_clone->Draw("same l p");

  TString outFileName = "outkmaxalfa.root";
  TFile *outFile = new TFile(outFileName, "recreate");
  outFile->cd();
  cout << endl << "Output in " << outFileName.Data() << endl << endl;
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();
  c7->Write();
  outFile->Close();
}
