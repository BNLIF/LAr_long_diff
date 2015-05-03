void plot_stable(){
  ifstream infile("./results/result.txt");
  Double_t range_min[200],range_max[200];
  Double_t fit_diff[200], fit_diff_err[200];
  Double_t chi2ndf[200],vel[200],vel_err[200];
  Int_t temp;
  for (Int_t i=0;i!=185;i++){
    infile >> temp;
    infile >> range_min[temp] >> range_max[temp] 
	   >> fit_diff[temp] >> fit_diff_err[temp]
	   >> chi2ndf[temp] >> vel[temp] >> vel_err[temp];
    
    vel_err[temp]*=sqrt(chi2ndf[temp]);
    fit_diff_err[temp]*=sqrt(chi2ndf[temp]);
  }
  
  Double_t x[17];
  for (Int_t i=0;i!=17;i++){
    x[i] = i+1;
  }
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  vel[185] = vel[52];
  vel_err[185] = vel_err[52];
  TGraphErrors *g1 = new TGraphErrors(17,x,&vel[169],0,&vel_err[169]);
  g1->Draw("A*");
  g1->GetXaxis()->SetTitle("Setting");
  g1->SetTitle("Drift Time (ns)");
  
  cout << vel[52] << " " << vel_err[52] << endl;

  c1->cd(2);
  fit_diff[185] = fit_diff[52];
  fit_diff_err[185] = fit_diff_err[52];
  TGraphErrors *g1 = new TGraphErrors(17,x,&fit_diff[169],0,&fit_diff_err[169]);
  g1->Draw("A*");
  g1->GetXaxis()->SetTitle("Setting");
  g1->SetTitle("Diffusion Time (ns)");
}
