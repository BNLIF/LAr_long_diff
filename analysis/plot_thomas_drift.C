void plot_thomas_drift(){
  TFile *file = new TFile("result_thomas_new.root");
  TTree *T = (TTree*)file->Get("T");
  Int_t fileno;
  Double_t E, dis, diff_time, drift_time;
  Double_t cutoff, m_charge;
  Double_t m_diff_time, m_diff_time_err;
  Double_t m_drift_time;

  T->SetBranchAddress("fileno",&fileno);
  T->SetBranchAddress("E",&E);
  T->SetBranchAddress("dis",&dis);
  T->SetBranchAddress("diff_time",&diff_time);
  T->SetBranchAddress("drift_time",&drift_time);
  T->SetBranchAddress("cutoff",&cutoff);
  T->SetBranchAddress("m_charge",&m_charge);
  T->SetBranchAddress("m_diff_time",&m_diff_time);
  T->SetBranchAddress("m_diff_time_err",&m_diff_time_err);
  T->SetBranchAddress("m_drift_time",&m_drift_time);
  
  Int_t ncount[8]={0,0,0,0,0,0,0,0};
  Double_t norm_drift_time[8][100],mE[8][100];
  Double_t m_norm_drift_time[8][100];

  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (dis!=0){
      Int_t n = 0;
      if (dis==7.5) n=1;
      if (dis==10) n=2;
      if (dis==15) n=3;
      if (dis==20) n=4;
      if (dis==30) n=5;
      if (dis==45) n=6;
      if (dis==60) n=7;
      
      
      norm_drift_time[n][ncount[n]] = drift_time/dis;
      m_norm_drift_time[n][ncount[n]] = m_drift_time/dis;
      mE[n][ncount[n]] = E;
      ncount[n]++;
    }
  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(3,3);
  c1->cd(1);
  c1_1->SetLogx(1);
  TGraph *g1 = new TGraph(ncount[0],mE[0],m_norm_drift_time[0]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("5 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[0],mE[0],norm_drift_time[0]);
  g2->Draw("Lsame");
  g2->SetLineColor(4);

  c1->cd(2);
  Int_t n = 1;
  TGraph *g1 = new TGraph(ncount[1],mE[1],m_norm_drift_time[1]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("7.5 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[1],mE[1],norm_drift_time[1]);
  g2->Draw("Lsame");
  g2->SetLineColor(4);

  c1->cd(3);
  Int_t n = 2;
  TGraph *g1 = new TGraph(ncount[n],mE[n],m_norm_drift_time[n]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("10 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[n],mE[n],norm_drift_time[n]);
  g2->Draw("Lsame");
g2->SetLineColor(4);

  c1->cd(4);
  Int_t n = 3;
  TGraph *g1 = new TGraph(ncount[n],mE[n],m_norm_drift_time[n]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("15 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[n],mE[n],norm_drift_time[n]);
  g2->Draw("Lsame");
g2->SetLineColor(4);
  c1->cd(5);
  Int_t n = 4;
  TGraph *g1 = new TGraph(ncount[n],mE[n],m_norm_drift_time[n]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("20 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[n],mE[n],norm_drift_time[n]);
  g2->Draw("Lsame");
g2->SetLineColor(4);
  c1->cd(6);
  Int_t n = 5;
  TGraph *g1 = new TGraph(ncount[n],mE[n],m_norm_drift_time[n]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("30 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[n],mE[n],norm_drift_time[n]);
  g2->Draw("Lsame");
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(4);
  g2->SetLineColor(4);

  c1->cd(7);
  Int_t n = 6;
  TGraph *g1 = new TGraph(ncount[n],mE[n],m_norm_drift_time[n]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("45 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[n],mE[n],norm_drift_time[n]);
  g2->Draw("Lsame");
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(4);
  g2->SetLineColor(4);

  c1->cd(8);
  Int_t n = 7;
  TGraph *g1 = new TGraph(ncount[n],mE[n],m_norm_drift_time[n]);
  g1->Draw("A*");
  g1->SetMarkerStyle(20);
  g1->SetTitle("60 mm");
  g1->GetYaxis()->SetTitle("1/v_{Drift}");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  TGraph *g2 = new TGraph(ncount[n],mE[n],norm_drift_time[n]);
  g2->Draw("*same");
  g2->SetMarkerStyle(21);
  g2->SetMarkerColor(4);
  g2->SetLineColor(4);

  

}
