void slices(int posflag, int verflag)
{
  // make slices from T2res as a function of energy
 
  TFile *f1;
  if( verflag ==0) {
    f1 = new TFile("prova_sel100.root");
  } else if( verflag ==1){
    f1 = new TFile("prova_sel200.root");
  } else if( verflag == 2){
    f1 = new TFile("prova.root");
  }else{
    cout << " VerFlag = " << verflag << " Not Yet Implemented !!!" << endl;
    exit(0);
  }
  
  f1->ls();

  TCanvas *c0;
  c0= new TCanvas("c0","Mips vs Z");
  c0->Divide(3,3);
  TH2F *EZ_pos = (TH2F*)f1->Get("EvsZ_mpos");
  TH2F *EZ_neg = (TH2F*)f1->Get("EvsZ_mneg");
  TH1F *Emip_pos[10];
  TH1F *Emip_neg[10];
		

  f1->cd();
  for(int j=0; j<9; j++){
    int bin_min=10+2*j;
    int bin_max=10+2*(j+1);
    c0->cd(j+1);
    if( posflag==0){
      Emip_pos[j] = (TH1F*)EZ_pos->ProjectionY(Form("Emip_pos%d",j),bin_min,bin_max,"");
      Emip_pos[j]->Draw("e");
    } else if ( posflag ==1){
      Emip_neg[j] = (TH1F*)EZ_neg->ProjectionY(Form("Emip_neg%d",j),bin_min,bin_max,"");     
      Emip_neg[j]->Draw("e");
    }

  }
   
  return;

}
