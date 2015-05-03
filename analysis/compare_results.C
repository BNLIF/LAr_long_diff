void compare_results(){
  ifstream infile1("./results/result.txt");
  ifstream infile2("./results/result_final.txt");

  Double_t r_fileno[200];
  Double_t r_fit_min[200],r_fit_max[200];
  Double_t r_fit_rms[200],r_fit_mean[200];
  Double_t r_fit_rms_err[200],r_fit_mean_err[200];
  Double_t temp,chi2;
  Double_t temp1,temp2;
  Int_t fit_flag=0;
  for (Int_t i=0;i!=184;i++){
    infile1 >> r_fileno[i];
    infile2 >> temp;
    Int_t n = r_fileno[i];
    if (n<200 && n>=0){
      infile1 >> r_fit_min[n] >> r_fit_max[n];
      infile1 >> r_fit_rms[n] >> r_fit_rms_err[n] >> chi2;
      infile1 >> r_fit_mean[n] >> r_fit_mean_err[n];
      
      r_fit_rms_err[n] = fabs(r_fit_rms_err[n]*sqrt(chi2));
      r_fit_mean_err[n] = fabs(r_fit_mean_err[n]*sqrt(chi2));

      infile2 >> temp >> temp;
      infile2 >> temp1 >> temp >> temp;
      infile2 >> temp2 >> temp;
    
      if (temp1 !=0){
	r_fit_rms[n] = (fabs(r_fit_rms[n]) - fabs(temp1))/temp1;
	r_fit_rms_err[n] /= temp1;
	r_fit_mean[n] = (r_fit_mean[n] - temp2)/temp2;
	r_fit_mean_err[n] /= temp2;

	if (fabs(r_fit_rms[n]) > r_fit_rms_err[n]) 
	  cout << n  << " " << temp1 << " " <<r_fit_rms[n]*temp1 + temp1<< endl;
      }else{
	r_fit_rms[n] =0;
	r_fit_rms_err[n] = 0.;
      }
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",800,400);
  c1->Divide(2,1);
  c1->cd(1);
  TGraphErrors *g1 = new TGraphErrors(185,r_fileno,r_fit_rms,0,r_fit_rms_err);
  g1->Draw("A*");
  c1->cd(2);
  TGraphErrors *g2 = new TGraphErrors(185,r_fileno,r_fit_mean,0,r_fit_mean_err);
  g2->Draw("A*");
}
