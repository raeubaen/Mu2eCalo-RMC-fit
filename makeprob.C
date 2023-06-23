void makeprob(int flag, float x0)
{
  TFile *f1 = new TFile("Reso.root");
  f1->ls();
  
  //f1->cd("ALL_cuts");
  //f1->ls();
  TH1F* hh[20];
  double_t ProbReso[20][300];
  float pgam;
  int jselgam;
  
  for (int jgam=0; jgam<20; jgam++){
    hh[jgam] = (TH1F*)f1->Get(Form("ALL_cuts/Reso_cutALL_%d",jgam));
    cout << " --- jgam --- " << jgam << endl;    
    //int jselgam = pgam/5;
    if( flag==1) {
      jselgam  = jgam;
      pgam = jgam*5;
    }
    if( jgam == jselgam){
      //hh[jgam]->Draw();
      double NevTot = hh[jgam]->Integral();
      int NumBin = hh[jgam]->GetNbinsX();
      float xmin = hh[jgam]->GetXaxis()->GetXmin();
      float xmax = hh[jgam]->GetXaxis()->GetXmax();
      float step = (xmax-xmin)/NumBin;
      //cout << " Step " << step << endl;
      
      double xval[300];
      double yval[300];
      double sum=0;
      double distance = x0-pgam;
      int xbintobefound = (distance-xmin)/step;
      // cout << "Xmin-Xmax-Step " << xmin << " " << xmax << " " << step << " X0-Pgam " << distance << " IBIN  " << xbintobefound << endl;
      
      for(int j=0;j<NumBin;j++){
	if( NevTot ==0){
	  ProbReso[jgam][j] = 0;
	}else{
	  xval[j] = hh[jgam]->GetBinCenter(j+1);
	  yval[j] = hh[jgam]->GetBinContent(j+1);
	  ProbReso[jgam][j]= yval[j]/NevTot;
	}
	//sum = sum + yval[j];
      }
      //cout << " sum " << sum << " NevTot " << NevTot << endl;

      cout <<  " X0 " << x0 << " pgam " << pgam << " X0-pg " << distance << " IBIN  " << xbintobefound << " Prob " << ProbReso[jgam][xbintobefound] <<  endl;
      
    }
  }
  return;
}
