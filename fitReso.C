{
  TFile *inFile = new TFile("Reso.root");
  int firstslice = 6;
  int lastslice = 19;
  int nslice = (lastslice - firstslice);

  double *mpv = new double[nslice];
  double *mpv_err = new double[nslice];
  double *scale = new double[nslice];
  double *scale_err = new double[nslice];
  double *x = new double[nslice];

  TFile *outFile = new TFile("outfitReso.root", "recreate");
  for(int iSl=firstslice; iSl<lastslice; iSl++){
    TH1F *h = (TH1F*)inFile->Get(Form("ALL_cuts/Reso_cutALL_%i", iSl));
    TF1 f("f", "[0]*TMath::Landau(-x, [1], [2])", -40, 10);
    f.SetParameters(400, 5, 2);
    h->Fit(&f, "R");
    int iArr = iSl - firstslice;
    mpv[iArr] = f.GetParameter(1);
    mpv_err[iArr] = f.GetParError(1);
    scale[iArr] = f.GetParameter(2);
    scale_err[iArr] = f.GetParError(2);
    outFile->cd();
    h->Write();
    x[iArr] = iSl*5 + 2.5;
    cout << iSl << endl;
  }

  TGraphErrors *g_mpv = new TGraphErrors(nslice, x, mpv, 0, mpv_err);
  TGraphErrors *g_scale = new TGraphErrors(nslice, x, scale, 0, scale_err);
  g_mpv->SetName("MPV");
  g_scale->SetName("Scale");
  g_mpv->Fit("pol1");
  g_mpv->SetMarkerStyle(22);
  g_scale->Fit("pol1");
  g_scale->SetMarkerStyle(22);
  g_mpv->Write();
  g_scale->Write();

  double *mpv_par = g_mpv->GetFunction("pol1")->GetParameters();
  double *scale_par = g_scale->GetFunction("pol1")->GetParameters();
  cout << mpv_par[0] << " " << mpv_par[1] << endl;
  cout << scale_par[0] << " " << scale_par[1] << endl;

  double ethmax = 150;
  double slicew = 0.5;
  int nsliceout = (int)(ethmax/slicew);
  TF1 *mpv_f = new TF1("mpv_f", "pol1", -100, 100);
  TF1 *scale_f = new TF1("scale_f", "pol1", -100, 100);
  mpv_f->SetParameters(mpv_par);
  scale_f->SetParameters(scale_par);
  double fl_lim[2] = {-60, 20};
  TF1 *fluct = new TF1(
    "fl", "[0]*TMath::Landau(-x, [1], [2])",
    fl_lim[0], fl_lim[1]
  );

  outFile->mkdir("ALL_cuts");
  outFile->cd("ALL_cuts");
  double binw = 0.5;
  int bins = (int)((fl_lim[1]-fl_lim[0])/binw);
  TH1F *hist;
  for(int i=0; i<nsliceout; i++){
    double en = slicew/2 + i*slicew;
    double par[3] = {400, mpv_f->Eval(en), scale_f->Eval(en)};
    fluct->SetParameters(par); //Form("Reso_cutALL_%i", i),
    hist = new TH1F(Form("Reso_cutALL_%i", i), "fluct", bins, fl_lim[0], fl_lim[1]);
    hist->FillRandom("fl", 1e4);
    hist->Write();
  }
}
