

void plot_diff_all(){
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

    fit_diff_err[temp]*=sqrt(chi2ndf[temp]);
  }

  Int_t ncount[18]={0,0,0,0,0,0,0,0,0,
		    0,0,0,0,0,0,0,0,0};
  Double_t mE[18][100],mdis[18];
  Double_t m_norm_diff_time[18][100];
  Double_t m_norm_diff_time_err[18][100];
  Double_t m_norm_drift_time[18][100];
  Double_t m_norm_drift_time_err[18][100];

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
          
      //GAr
      mdis[n] = dis;
      mE[n][ncount[n]] = E ;
      
      m_norm_drift_time[n][ncount[n]] = fit_diff[fileno];
      m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(fit_diff_err[fileno],2) + pow(fit_diff[fileno]*0.05,2));


      ncount[n]++;
    }else if (i>=158 && i<169){
      if (dis==15){
	n=7;
      }else if (dis==30){
	n=8;
      }
     

      //SAr
      mdis[n] = dis;
      mE[n][ncount[n]] = E;
      m_norm_drift_time[n][ncount[n]] = fit_diff[fileno];
      m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(fit_diff_err[fileno],2) + pow(fit_diff[fileno]*0.035,2));
      
      if (vel_err[fileno] >0)
	ncount[n]++;
    }else if (i<=97&&i<158){
      if (dis==5){
	n=9;
      }else if (dis==7.5){
	n=10;
      }else if (dis==10){
	n=11;
      }else if (dis==15){
	n=12;
      }else if (dis==20){
      	n=13;
      }else if (dis==30){
	n=14;
      }else if (dis==45){
	n=15;
      }else if (dis==60){
	n=16;
      }
      

      //LAr
      mdis[n] = dis;
      mE[n][ncount[n]] = E;
      
     Double_t dE = 0.025*sqrt(5./dis);

     //cout << argon.vDrift(temperature,E)/argon.vDrift(temperature,E+dE) << endl;
     
     m_norm_drift_time[n][ncount[n]] = fit_diff[fileno];
     m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(fit_diff_err[fileno],2) + pow(fit_diff[fileno]*0.035,2));

     
     if (vel_err[fileno] >0)
      ncount[n]++;
    }
  }


  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->SetFillColor(10);
  c1->Divide(3,1);
  c1->cd(1);
  c1_1->SetLogx(1);
  c1_1->SetLogy(1);

  
  TH2F *h1 = new TH2F("h1","h1",100,0.005,3.5,100,50.,3000);
  h1->Draw();
  TGraphErrors **g1 = new TGraphErrors*[16];
  TLegend *le1 = new TLegend(0.6,0.3,0.89,0.89);
  le1->SetFillColor(10);
  for (Int_t i=0;i!=7;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_drift_time[i],0,m_norm_drift_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i);
    g1[i]->SetMarkerColor(1+i);
    if (i==4){
      g1[i]->SetMarkerColor(8);
    }

    g1[i]->SetMarkerSize(1.2);
    if (i<7){
      le1->AddEntry(g1[i],Form("%2.0f mm",mdis[i]),"p");
    }else{
      le1->AddEntry(g1[i],Form("%2.0f mm",mdis[i]),"p");
    }
  }
 
  le1->Draw();
  h1->SetXTitle("E (V/cm)");
  h1->SetYTitle("T_{Diffusion} (ns)");
  h1->SetTitle("GAr");
  h1->GetYaxis()->SetTitleOffset(1.2);
   
  c1->cd(2);
  c1_2->SetLogx(1);
  c1_2->SetLogy(1);

  TH2F *h2 = new TH2F("h2","h2",100,0.005,4.5,100,40.,10000);
  h2->Draw();
  TLegend *le2 = new TLegend(0.6,0.3,0.89,0.89);
  le2->SetFillColor(10);
  for (Int_t i=9;i!=17;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_drift_time[i],0,m_norm_drift_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i-9);
    g1[i]->SetMarkerColor(1+i-9);
    g1[i]->SetMarkerSize(1.2);
    if (i==13)
      g1[i]->SetMarkerColor(4);
    le2->AddEntry(g1[i],Form("%2.1f mm",mdis[i]),"p");
  }
 
  le2->Draw();
  h2->SetXTitle("E (kV/cm)");
  h2->SetYTitle("");
  h2->SetTitle("LAr");



  c1->cd(3);

  TH2F *h3 = new TH2F("h3","h3",100,0.005,1.3,100,0.,1600);
  h3->Draw();
  TLegend *le3 = new TLegend(0.6,0.3,0.89,0.89);
  le3->SetFillColor(10);
  for (Int_t i=7;i!=9;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_drift_time[i],0,m_norm_drift_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i-7);
    g1[i]->SetMarkerColor(1+i-7);
    g1[i]->SetMarkerSize(1.2);
    if (i<7){
      le3->AddEntry(g1[i],Form("%2.0f mm",mdis[i]),"p");
    }else{
      le3->AddEntry(g1[i],Form("%2.0f mm",mdis[i]),"p");
    }
  }
  le3->Draw();
  h3->SetXTitle("E (kV/cm)");
  h3->SetYTitle("");
  h3->SetTitle("SAr");


 
}
