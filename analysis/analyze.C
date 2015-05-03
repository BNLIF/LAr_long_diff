void analyze(){
  TFile *file = new TFile("merge.root");
  
  TTree *T = (TTree*)file->Get("T");
  Int_t fileno;
  Double_t E, dis, diff_time, drift_time;
  Double_t cutoff, m_charge;
  Double_t m_diff_time, m_diff_time_err;
  Double_t m_drift_time,m_drift_time_err;
  
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
  T->SetBranchAddress("m_drift_time_err",&m_drift_time_err);

  Double_t x2ndf;
  Double_t min,max;
  T->SetBranchAddress("x2ndf",&x2ndf);
  T->SetBranchAddress("min",&min);
  T->SetBranchAddress("max",&max);

  TH1F *h1 = new TH1F("h1","h1",1000,2,0.6);
  Int_t ncount[100];
  Double_t Efield[100];
  for (Int_t i=0;i!=100;i++){
    ncount[i] = 0;
  }

  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (E>0&&dis!=45 ){
      h1->Fill(log10(E));
    }
  }
  h1->Draw();

  Int_t n=0;
  for (Int_t i=0;i!=h1->GetNbinsX();i++){
    ncount[n] = h1->GetBinContent(i+1);
    Efield[n] = pow(10,h1->GetBinCenter(i+1));
    if (ncount[n]>=1){
      n++;
    }
  }
  cout << n << endl;
  
  TGraphErrors **g1 = new TGraphErrors*[n];
  for (Int_t i=0;i!=n;i++){
    g1[i] = new TGraphErrors();
  }

  for (Int_t i=0;i!=100;i++){
    ncount[i] = 0;
  }

  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (E>0 && dis!=45&&dis!=0&&m_drift_time!=0&&E>=0.1){
      for (Int_t j=0;j!=n;j++){
	//cout << log10(E) << " " << log10(Efield[j]) << " " << dis  << endl;
	if (fabs(log10(E)-log10(Efield[j]))<3./1000.){
	  cout << j << " " << Efield[j] << " " << dis << " " << 
	    m_diff_time << " " << m_diff_time_err*sqrt(x2ndf) << 
	    " " << m_drift_time << " " << m_drift_time_err << endl;
	  
	  // g1[j]->SetPoint(ncount[j],dis,m_diff_time*m_diff_time);
	  // g1[j]->SetPointError(ncount[j],0,2*m_diff_time*m_diff_time_err*sqrt(x2ndf));
	  // ncount[j]++;
	} 
      }
    }
  }
  
  // TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  // c1->Divide(4,4);
  // for (Int_t i=0;i!=n;i++){
  //   c1->cd(i+1);
  //   // cout << g1[i]->GetN() << endl;
  //   g1[i]->Draw("A*");
  //   g1[i]->Fit("pol1");
  //   g1[i]->SetTitle(Form("E-field %3.2f",Efield[i]));
  // }
  

}
