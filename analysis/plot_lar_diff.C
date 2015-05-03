void plot_lar_diff(){
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
  Double_t mE[9][100],mdis[9];
  Double_t m_norm_diff_time[9][100];
  Double_t m_norm_diff_time_err[9][100];
  Double_t m_norm_drift_time[9][100];
  Double_t m_norm_drift_time_err[9][100];

  for (Int_t i=0;i!=t1->GetEntries();i++){
    t1->GetEntry(i);
    int n;
    if (i<=97&&i<158){
      if (dis==5){
	n=0;
      }else if (dis==7.5){
	n=1;
      }else if (dis==10){
	n=2;
      }else if (dis==15){
	n=3;
      }else if (dis==20){
      	n=4;
      }else if (dis==30){
	n=5;
      }else if (dis==45){
	n=6;
      }else if (dis==60){
	n=7;
      }
      
      mdis[n] = dis;
      mE[n][ncount[n]] = E;
      m_norm_drift_time[n][ncount[n]] = pow(fit_diff[fileno]/vel[fileno],2)/2.*dis/10.*E*1000.;
      m_norm_drift_time_err[n][ncount[n]] = m_norm_drift_time[n][ncount[n]]*2*sqrt(pow(vel_err[fileno]/vel[fileno],2)+pow(fit_diff_err[fileno]/fit_diff[fileno],2));
      ncount[n]++;
      //      cout << dis << " " << E << " " << fileno << endl;
    }
  }


  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->SetLogx(1);
  c1->SetLogy(1);
  TH2F *h1 = new TH2F("h1","h1",100,0.001,3.5,100,0.001,10.0);
  h1->Draw();
  TGraphErrors **g1 = new TGraphErrors*[9];
  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  for (Int_t i=0;i!=8;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_drift_time[i],0,m_norm_drift_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i);
    g1[i]->SetMarkerColor(1+i);
    g1[i]->SetMarkerSize(1.5);
    le1->AddEntry(g1[i],Form("LAr %2.1f mm",mdis[i]),"p");

  }
  le1->Draw();
  h1->SetXTitle("E (kV/cm)");
  h1->SetYTitle("eV");
  h1->SetTitle("Electron Energy LAr");
  
  
}
