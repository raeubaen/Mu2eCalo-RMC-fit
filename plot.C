{
    TFile *f0 = new TFile("kmax90alpha2.root");
    TFile *f1 = new TFile("kmax90.2alpha1.9.root");

    TCanvas *c0 = (TCanvas*)f0->Get("c4");
    TCanvas *c1 = (TCanvas*)f1->Get("c4");

    TH1F *h0 = (TH1F*)(c0->GetListOfPrimitives()->FindObject("eexp"));
    TH1F *h0r = (TH1F*)(c0->GetListOfPrimitives()->FindObject("ereco"));
    TH1F *h1 = (TH1F*)(c1->GetListOfPrimitives()->FindObject("eexp"));
    TH1F *h1r = (TH1F*)(c1->GetListOfPrimitives()->FindObject("ereco"));

    TCanvas *c = new TCanvas("c", "c");
    h0->SetLineColor(kBlack);
    h0->Draw();
    h0r->SetLineColor(kRed);
    h0r->Draw("same");
    h1->SetLineColor(kBlue);
    h1->Draw("same");
    h1r->SetLineColor(kViolet);
    h1r->Draw("same");

    TLegend *l = new TLegend();
    l->AddEntry(h0, "Conv. - kmax=90, alpha=2");
    l->AddEntry(h0r, "MC. - kmax=90, alpha=2");
    l->AddEntry(h1, "Conv. - kmax=90.2, alpha=1.9");
    l->AddEntry(h1r, "MR. - kmax=90.2, alpha=1.9");

    l->Draw();
}