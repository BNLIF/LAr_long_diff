void analyze_update(Int_t set = 0){

  Int_t setting[200];
  Double_t E[200], dis[200];
  Double_t m_diff_time[200], m_diff_time_err[200];
  Double_t m_drift_time[200],m_drift_time_err[200];

  ifstream infile("result_update.dat");
  for (Int_t i=0;i!=163;i++){
    infile >> setting[i] >> E[i] >> dis[i] 
	   >> m_diff_time[i] >> m_diff_time_err[i]
	   >> m_drift_time[i] >> m_drift_time_err[i];
  }
  
  TH1F *h1 = new TH1F("h1","h1",1000,-1.1,0.65);
  Int_t ncount[100];
  Double_t Efield[100];
  for (Int_t i=0;i!=100;i++){
    ncount[i] = 0;
  }

  for (Int_t i=0;i!=163;i++){
    if (setting[i] == set){
      if (E[i]>0.099){
	h1->Fill(log10(E[i]));
      }
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
  //cout << n << endl;
  
 

  for (Int_t i=0;i!=100;i++){
    ncount[i] = 0;
  }

  ofstream outfile;
  if (set==0){
    outfile.open("fit_result_gar.dat");
  }else if (set==1){
    outfile.open("fit_result_lar.dat");
  }else if (set==2){
    outfile.open("fit_result_sar.dat");
  }

  for (Int_t i=0;i!=163;i++){
    if (setting[i] == set && E>0.099 ){
      for (Int_t j=0;j!=n;j++){
       	if (fabs(log10(E[i])-log10(Efield[j]))<3./1000.){
   	  outfile << j << " " << Efield[j] << " " << dis[i] << " " << 
   	    m_diff_time[i] << " " << m_diff_time_err[i] << 
   	    " " << m_drift_time[i] << " " << m_drift_time_err[i] << endl;
	  

  // 	  // ncount[j]++;
  	} 
      }
    }
  }
  
  // // TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  // // c1->Divide(4,4);
  // // for (Int_t i=0;i!=n;i++){
  // //   c1->cd(i+1);
  // //   // cout << g1[i]->GetN() << endl;
  // //   g1[i]->Draw("A*");
  // //   g1[i]->Fit("pol1");
  // //   g1[i]->SetTitle(Form("E-field %3.2f",Efield[i]));
  // // }
  

}
