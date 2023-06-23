{
    TFile *inFile = new TFile("Reso.root");
    TFile *outFile = new TFile("gaussReso.root", "recreate");

    outFile->mkdir("ALL_cuts");
    outFile->cd("ALL_cuts");
    for(int i=0; i<19; i++){
        TH1F *h = new TH1F(Form("Reso_cutALL_%i", i), "h", 300, -10, 10);
        h->FillRandom("gaus", 1e6);
        h->Write();
    }
}
