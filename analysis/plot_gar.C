void plot_gar(){
  Double_t dis, E;
  Int_t fileno;
  
  TChain *t1 = new TChain("T","T");
  t1->AddFile("thomas.root");
  t1->SetBranchAddress("dis",&dis);
  t1->SetBranchAddress("E",&E);
  t1->SetBranchAddress("fileno",&fileno);

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

  Int_t ncount[9]={0,0,0,0,0,0,0,0,0};
  Double_t mE[9][100];
  Double_t m_norm_diff_time[9][100];
  Double_t m_norm_diff_time_err[9][100];
  Double_t m_norm_drift_time[9][100];
  Double_t m_norm_drift_time_err[9][100];

  for (Int_t i=0;i!=t1->GetEntries();i++){
    t1->GetEntry(i);
    int n;
    if (i>=100&&i<158){
      if (dis==5){
	n=0;
      }else if (dis==10){
	n=1;
      }else if (dis==15){
	n=2;
      }else if (dis==20){
	n=3;
      }else if (dis==30){
      	n=4;
      }else if (dis==45){
	n=5;
      }else if (dis==60){
	n=6;
      }
      
      mE[n][ncount[n]] = E;
      m_norm_drift_time[n][ncount[n]] = vel[fileno];
      m_norm_drift_time_err[n][ncount[n]] = vel_err[fileno];
      m_norm_diff_time[n][ncount[n]] = fit_diff[fileno];
      m_norm_diff_time_err[n][ncount[n]] = fit_diff_err[fileno];
      ncount[n]++;
      //      cout << dis << " " << E << " " << fileno << endl;
    }else if (i>=158 && i<169){
      if (dis==15){
	n=7;
      }else if (dis==30){
	n=8;
      }
      mE[n][ncount[n]] = E;
      m_norm_drift_time[n][ncount[n]] = vel[fileno];
      m_norm_drift_time_err[n][ncount[n]] = vel_err[fileno];
      m_norm_diff_time[n][ncount[n]] = fit_diff[fileno];
      m_norm_diff_time_err[n][ncount[n]] = fit_diff_err[fileno];
      if (vel_err[fileno] >0)
	ncount[n]++;
    }
  }


  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(3,3);
  for (Int_t i=0;i!=9;i++){
    c1->cd(i+1);
    TGraphErrors *g1 = new TGraphErrors(ncount[i],mE[i],m_norm_drift_time[i],0,m_norm_drift_time_err[i]);
    g1->Draw("A*");
  }

  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  c2->Divide(3,3);
  for (Int_t i=0;i!=9;i++){
    c2->cd(i+1);
    TGraphErrors *g1 = new TGraphErrors(ncount[i],mE[i],m_norm_diff_time[i],0,m_norm_diff_time_err[i]);
    g1->Draw("A*");
  }

}
