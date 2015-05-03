void plot_raw_thomas(){
  TFile *file = new TFile("thomas.root");
  TH1D *h5 = (TH1D*)file->Get("S1_45");
  TH1D *h10 = (TH1D*)file->Get("S1_2"); 
  TH1D *h15 = (TH1D*)file->Get("S1_17");
 
  TH1D *h30 = (TH1D*)file->Get("S1_23");
  TH1D *h45 = (TH1D*)file->Get("S1_35");
  TH1D *h60 = (TH1D*)file->Get("S1_57");

  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(3,2);
  c1->cd(1);
  h5->Draw();
  h5->SetTitle("5 mm 0.1 kV/cm");

  c1->cd(2);
  h10->Draw();
  h10->SetTitle("10 mm 0.2 kV/cm");

  c1->cd(3);
  h15->Draw();
  h15->SetTitle("15 mm 1.5 kV/cm");

  c1->cd(4);
  h30->Draw();
  h30->SetTitle("30 mm 0.2 kV/cm");

  c1->cd(5);
  h45->Draw();
  h45->SetTitle("45 mm 0.2 kV/cm");

  c1->cd(6);
  h60->Draw();
  h60->SetTitle("60 mm 0.5 kV/cm");
}
