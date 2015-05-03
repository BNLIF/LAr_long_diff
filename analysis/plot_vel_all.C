#include "LAr.C"

void plot_vel_all(){
  Double_t dis, E;
  Int_t fileno;
  
   LAr argon;
  Double_t temperature = 87;

  TChain *t1 = new TChain("T","T");
  t1->AddFile("thomas.root");
  t1->SetBranchAddress("dis",&dis);
  t1->SetBranchAddress("E",&E);
  t1->SetBranchAddress("fileno",&fileno);

  ifstream infile("./results/result_final.txt");
  ifstream infile1("./results/result.txt");
  Double_t range_min[200],range_max[200];
  Double_t fit_diff[200], fit_diff_err[200];
  Double_t chi2ndf[200],vel[200],vel_err[200];
  Int_t temp;
  Double_t temp1, temp2,temp3;
  for (Int_t i=0;i!=185;i++){
    infile >> temp;
    infile >> range_min[temp] >> range_max[temp] 
	   >> fit_diff[temp] >> fit_diff_err[temp]
	   >> chi2ndf[temp] >> vel[temp] >> vel_err[temp];
    
    infile1 >> temp3 >> temp3 >> temp3;
    infile1 >> temp1 >> temp3 >> temp3 >> temp2 >> temp3;

    //subtract time 0, same for everybody
    vel[temp] -= 27;
    vel_err[temp]*=sqrt(chi2ndf[temp]);
    fit_diff_err[temp]*=sqrt(chi2ndf[temp]);

    vel_err[temp] = sqrt(pow(vel_err[temp],2) + 4 + pow(vel[temp] - temp2,2));
    fit_diff_err[temp] = sqrt(pow(fit_diff_err[temp],2) + pow(fit_diff[temp]-temp1,2));
  }

  Int_t ncount[18]={0,0,0,0,0,0,0,0,0,
		    0,0,0,0,0,0,0,0,0};
  Double_t mE[18][100],mdis[18];
  Double_t m_norm_diff_time[18][100];
  Double_t m_norm_diff_time_err[18][100];
  Double_t m_norm_drift_time[18][100];
  Double_t m_norm_drift_time_err[18][100];

  ofstream outfile("result_update.dat");

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
      
      Double_t dE = 0.025*sqrt(5./dis);
      Double_t p0 = 0.34953;
      Double_t p1 = 0.07498;
      Double_t p2 = 6.6984e-3;
      //GAr
      mdis[n] = dis;
      mE[n][ncount[n]] = E *1000./1134.8;
      m_norm_drift_time[n][ncount[n]] = dis/10./(vel[fileno]*1e-9-12e-9)/(E*1000.);
      m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(m_norm_drift_time[n][ncount[n]]*sqrt(pow(vel_err[fileno],2)+144)/vel[fileno],2) //stat
						 + pow(0.1/10./(vel[fileno]*1e-9)/(E*1000.),2) // distance
						 +pow(0.07*m_norm_drift_time[n][ncount[n]],2) //temperature
						 +pow(m_norm_drift_time[n][ncount[n]]*(1-(p0+p1*log(E)+p2*log(E)*log(E))/(p0+p1*log(E+dE)+p2*log(E+dE)*log(E+dE))),2));
      m_norm_diff_time[n][ncount[n]] = dis/10./(vel[fileno]*1e-9-12e-9)/(E*1000.)*E/1000.;
      m_norm_diff_time_err[n][ncount[n]] = sqrt(pow(m_norm_diff_time[n][ncount[n]]*sqrt(pow(vel_err[fileno],2)+144)/vel[fileno],2) //stat
						+ pow(0.1/10./(vel[fileno]*1e-9)/(E*1000.)*E/1000.,2) //distance
						+pow(0.07*m_norm_diff_time[n][ncount[n]],2) //temperature
						+pow(m_norm_diff_time[n][ncount[n]]*(1-(p0+p1*log(E)+p2*log(E)*log(E))/(p0+p1*log(E+dE)+p2*log(E+dE)*log(E+dE))),2));

      outfile << 0 << " " << E << " " << dis << " " << vel[fileno] - 12 << " "<< sqrt(pow(vel_err[fileno],2)+144 + pow(0.07*(vel[fileno] - 12),2)) << " " << fit_diff[fileno] << " " << sqrt(pow(fit_diff_err[fileno],2) + pow(fit_diff[fileno]*0.05,2)) << endl;
      

      ncount[n]++;
    }else if (i>=158 && i<169){
      if (dis==15){
	n=7;
      }else if (dis==30){
	n=8;
      }
       Double_t dE = 0.025*sqrt(5./dis);
      Double_t p0 = 0.1953;
      Double_t p1 = 0.5735;
      Double_t p2 = -0.2546;

      //SAr
      mdis[n] = dis;
      mE[n][ncount[n]] = E;
      m_norm_drift_time[n][ncount[n]] = dis/10./(vel[fileno]*1e-9-20e-9)/(E*1000.);
      m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(m_norm_drift_time[n][ncount[n]]*sqrt(pow(vel_err[fileno],2)+400)/vel[fileno],2) //stat + grid-anode corr.
						 + pow(0.1/10./(vel[fileno]*1e-9)/(E*1000.),2) //distance
						 +pow(0.05*m_norm_drift_time[n][ncount[n]],2)
						 +pow(m_norm_drift_time[n][ncount[n]]*(1-(p0+p1*E+p2*E*E)/(p0+p1*(E+dE)+p2*(E+dE)*(E+dE))),2));//temperature

      m_norm_diff_time[n][ncount[n]] = dis/10./(vel[fileno]*1e-9-20e-9)/(E*1000.)*E/1000.;
      m_norm_diff_time_err[n][ncount[n]] = sqrt(pow(m_norm_diff_time[n][ncount[n]]*sqrt(pow(vel_err[fileno],2)+400)/vel[fileno],2)  //stat  + grid-anode corr.
						+ pow(0.1/10./(vel[fileno]*1e-9)/(E*1000.)*E/1000.,2) //distance
						+pow(0.05*m_norm_diff_time[n][ncount[n]],2) // temperature
						+pow(m_norm_diff_time[n][ncount[n]]*(1-(p0+p1*E+p2*E*E)/(p0+p1*(E+dE)+p2*(E+dE)*(E+dE))),2));

      if (fit_diff_err[fileno]>0)
       outfile << 2 << " " << E << " " << dis << " " << vel[fileno] - 20 << " "<< sqrt(pow(vel_err[fileno],2)+400 + pow(0.05*(vel[fileno] - 20),2)) << " " << fit_diff[fileno] << " " << sqrt(pow(fit_diff_err[fileno],2) + pow(fit_diff[fileno]*0.035,2)) << endl;

      if (fit_diff_err[fileno] >0)
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
     
      m_norm_drift_time[n][ncount[n]] = dis/10./(vel[fileno]*1e-9 - 65e-9)/(E*1000.);
      m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(m_norm_drift_time[n][ncount[n]]*sqrt(pow(vel_err[fileno],2) + 50)/vel[fileno],2) //stat + drift vel uncertainties
						 +pow(0.1/10./(vel[fileno]*1e-9)/(E*1000.),2) //distance
						 +pow(0.05*m_norm_drift_time[n][ncount[n]],2) //temperature
						 +pow(m_norm_drift_time[n][ncount[n]]*(1-argon.vDrift(temperature,E)/argon.vDrift(temperature,E+dE)),2)
						 );

      m_norm_diff_time[n][ncount[n]] = dis/10./(vel[fileno]*1e-9 - 65e-9)/(E*1000.)*E/1000.;
      m_norm_diff_time_err[n][ncount[n]] = sqrt(pow(m_norm_diff_time[n][ncount[n]]*sqrt(pow(vel_err[fileno],2) + 50)/vel[fileno],2) //stat + drift vel uncertainties
						+ pow(0.1/10./(vel[fileno]*1e-9)/(E*1000.)*E/1000.,2) //distance
						+pow(0.05*m_norm_diff_time[n][ncount[n]],2) //temperature
						+pow(m_norm_diff_time[n][ncount[n]]*(1-argon.vDrift(temperature,E)/argon.vDrift(temperature,E+dE)),2));
      if (fit_diff_err[fileno]>0)
      outfile << 1 << " " << E << " " << dis << " " << vel[fileno] - 65 << " "<< sqrt(pow(vel_err[fileno],2)+50 + pow(0.05*(vel[fileno] - 65),2)) << " " << fit_diff[fileno] << " " << sqrt(pow(fit_diff_err[fileno],2) + pow(fit_diff[fileno]*0.035,2)) << endl;

      if (fit_diff_err[fileno] >0)
      ncount[n]++;
    }
  }


  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->SetFillColor(10);
  c1->Divide(3,2,0);
  c1->cd(1);
  c1_1->SetLogx(1);
  c1_1->SetLeftMargin(0.17);
  c1_4->SetLeftMargin(0.17);
  c1_4->SetBottomMargin(0.15);
  c1_1->SetRightMargin(0.01);
  c1_4->SetRightMargin(0.01);

  c1_2->SetLeftMargin(0.15);
  c1_5->SetLeftMargin(0.15);
  c1_5->SetBottomMargin(0.15);
  c1_2->SetRightMargin(0.01);
  c1_5->SetRightMargin(0.01);
  
  c1_3->SetLeftMargin(0.15);
  c1_6->SetLeftMargin(0.15);
  c1_6->SetBottomMargin(0.15);
  c1_3->SetRightMargin(0.01);
  c1_6->SetRightMargin(0.01);
  
  TH2F *h4 = new TH2F("h4","h4",100,0.005,3.5,100,0.05,0.65);
  h4->Draw();
  TGraphErrors **g1 = new TGraphErrors*[16];
  for (Int_t i=0;i!=7;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_diff_time[i],0,m_norm_diff_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i);
    g1[i]->SetMarkerColor(1+i);
    g1[i]->SetMarkerSize(1.2);
    if (i==4){
      g1[i]->SetMarkerColor(8);
    }
   
  }
  h4->SetXTitle("E (kV/cm)");
  h4->SetYTitle("v (cm/#mus)");
  h4->SetTitle("Gas Argon");
 
  
  c1->cd(2);
  c1_2->SetBottomMargin(0);
  c1_2->SetLogx(1);
  TH2F *h5 = new TH2F("h5","h5",100,0.005,4.5,100,-0.02,0.36);
  h5->Draw();
  for (Int_t i=9;i!=17;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_diff_time[i],0,m_norm_diff_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i-9);
    g1[i]->SetMarkerColor(1+i-9);
    g1[i]->SetMarkerSize(1.2);

    if (i==13)
      g1[i]->SetMarkerColor(4);
  }
  h5->SetXTitle("E (kV/cm)");
  h5->SetYTitle("");
  h5->SetTitle("Liquid Argon");



  c1->cd(3);
  c1_3->SetBottomMargin(0);
  TH2F *h6 = new TH2F("h6","h6",100,0.005,1.3,100,0.15,0.56);
  h6->Draw();
  for (Int_t i=7;i!=9;i++){
    g1[i] = new TGraphErrors(ncount[i],mE[i],m_norm_diff_time[i],0,m_norm_diff_time_err[i]);
    g1[i]->Draw("*same");
    g1[i]->SetMarkerStyle(20+i-7);
    g1[i]->SetMarkerColor(1+i-7);
    g1[i]->SetMarkerSize(1.2);
  }
  h6->SetXTitle("E (kV/cm)");
  h6->SetYTitle("");
  h6->SetTitle("Solid Argon");




  c1->cd(4);
  c1_4->SetLogx(1);
  c1_4->SetLogy(1);
  c1_4->SetTopMargin(0);
  c1_4->SetFrameBorderMode(0);
  TH2F *h1 = new TH2F("h1","h1",100,0.005,3.5,100,100.,17000);
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
  h1->SetXTitle("E/P (V/cm/torr)");
  h1->SetYTitle("#mu (cm^{2}/V/s)");
  h1->SetTitle("");
  h1->GetYaxis()->SetTitleOffset(1.2);
   
  c1->cd(5);
  c1_5->SetLogx(1);
  c1_5->SetTopMargin(0);
  TH2F *h2 = new TH2F("h2","h2",100,0.005,4.5,100,0.,2100);
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
  TLine *l1 = new TLine(0.005,518,0.05,518);
  l1->Draw();
  l1->SetLineWidth(1.5);
  le2->AddEntry(l1,"#mu(E=0 @ 87 K)","l");
  le2->Draw();
  h2->SetXTitle("E (kV/cm)");
  h2->SetYTitle("");
  h2->SetTitle("");



  c1->cd(6);
  c1_6->SetTopMargin(0);
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
  h3->SetTitle("");


 
}
