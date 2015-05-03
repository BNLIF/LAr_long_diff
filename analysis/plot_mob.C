#include "LAr.C"

void plot_mob(){
  TChain *T = new TChain("T","T");
  T->AddFile("merge.root");
  T->AddFile("result_lar6_run1.root");

  //  TTree *T = (TTree*)file->Get("T");
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
  
  Int_t ncount[9]={0,0,0,0,0,0,0,0,0};
  Double_t norm_drift_time[9][1000],mE[9][1000];
  
  Double_t m_norm_drift_time[9][1000];
  Double_t m_norm_drift_time_err[9][1000];

   LAr argon;
  Double_t temperature = 89;

  

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
      if (dis==6.5) n=8;
      if (E>0){
	
	Double_t dE = 0.025*sqrt(5./dis);
	m_drift_time = (m_drift_time-0.3/argon.vDrift(temperature,0.3/0.3*10.)*1000.);// / argon.vDrift(temperature,E)*argon.vDrift(temperature,E+dE) ;
	//cout << 0.3/argon.vDrift(temperature,0.3/0.3*10.)*1000. << " " << argon.vDrift(temperature,0.3/0.3*10.) << " " << argon.vDrift(temperature,0.5)<< endl;
	norm_drift_time[n][ncount[n]] = argon.vDrift(temperature,E)/argon.vDrift(temperature,E+dE);
	m_norm_drift_time[n][ncount[n]] = 1./(m_drift_time/(dis))/E*1e5;
	m_norm_drift_time_err[n][ncount[n]] = sqrt(pow(m_norm_drift_time[n][ncount[n]] * sqrt(pow(1-norm_drift_time[n][ncount[n]],2)+pow(40/m_drift_time,2)+0.01*0.01),2)) ;

	cout << fileno << " " << m_drift_time << " " << m_norm_drift_time_err[n][ncount[n]]/m_norm_drift_time[n][ncount[n]] << endl;

	mE[n][ncount[n]] = E;
	ncount[n]++;
      }
    }
  }
  Double_t x[100],y[3][100];
  for (Int_t i=0;i!=100;i++){
    x[i] = 4.0/100.*(i+0.5);
    y[0][i] = argon.vDrift(temperature,x[i])/x[i]*100.;
    y[1][i] = argon.vDrift(temperature-3,x[i])/x[i]*100.;
    y[2][i] = argon.vDrift(temperature+3,x[i])/x[i]*100.;
    //cout << y[i] << endl;
  }


  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  TH2F *h1 = new TH2F("h1","h1",100,0.001,4.0,100,0.001,2000);
  h1->Draw();
  h1->SetTitle("Electron Mobility");
  h1->SetXTitle("E field (kV/cm)");
  h1->SetYTitle("#mu (cm^{2}/V/s)");
  h1->GetYaxis()->SetTitleOffset(1.3);

  TGraph *g100 = new TGraph(100,x,y[0]);
  TGraph *g200 = new TGraph(100,x,y[1]);
  TGraph *g300 = new TGraph(100,x,y[2]);
  g100->Draw("Lsame");
  g100->SetLineColor(2);
  g100->SetLineWidth(2.);

  g200->Draw("Lsame");
  g200->SetLineColor(1);
  g200->SetLineWidth(2.);

  g300->Draw("Lsame");
  g300->SetLineColor(4);
  g300->SetLineWidth(2.);


  TGraphErrors *g10 = new TGraphErrors(ncount[0],mE[0],m_norm_drift_time[0],0,m_norm_drift_time_err[0]);
  TGraphErrors *g20 = new TGraphErrors(ncount[1],mE[1],m_norm_drift_time[1],0,m_norm_drift_time_err[1]);
  TGraphErrors *g30 = new TGraphErrors(ncount[2],mE[2],m_norm_drift_time[2],0,m_norm_drift_time_err[2]);
  TGraphErrors *g40 = new TGraphErrors(ncount[3],mE[3],m_norm_drift_time[3],0,m_norm_drift_time_err[3]);
  TGraphErrors *g50 = new TGraphErrors(ncount[4],mE[4],m_norm_drift_time[4],0,m_norm_drift_time_err[4]);
  TGraphErrors *g60 = new TGraphErrors(ncount[5],mE[5],m_norm_drift_time[5],0,m_norm_drift_time_err[5]);
  TGraphErrors *g70 = new TGraphErrors(ncount[6],mE[6],m_norm_drift_time[6],0,m_norm_drift_time_err[6]);
  TGraphErrors *g80 = new TGraphErrors(ncount[7],mE[7],m_norm_drift_time[7],0,m_norm_drift_time_err[7]);
  TGraphErrors *g90 = new TGraphErrors(ncount[8],mE[8],m_norm_drift_time[8],0,m_norm_drift_time_err[8]);

  g10->Draw("*same");
  g10->SetMarkerStyle(20);
  g10->SetMarkerColor(1); 
  g10->SetMarkerSize(1.0);
  g20->Draw("*same");
  g20->SetMarkerStyle(21);
  g20->SetMarkerColor(3);
  g30->Draw("*same");
  g30->SetMarkerStyle(22);
  g30->SetMarkerColor(4);
  g40->Draw("*same");
  g40->SetMarkerStyle(23);
  g40->SetMarkerColor(6);
  
  g50->Draw("*same");
  g50->SetMarkerStyle(24);
  g50->SetMarkerColor(1);

  g60->Draw("*same");
  g60->SetMarkerStyle(24);
  g60->SetMarkerColor(8);
  // g70->Draw("*same");
  // g70->SetMarkerStyle(25);
  // g70->SetMarkerColor(2);
  
  g80->Draw("*same");
  g80->SetMarkerStyle(25);
  g80->SetMarkerColor(2);

  g90->Draw("*same");
  g90->SetMarkerStyle(25);
  g90->SetMarkerColor(4);


 

  TLegend *le1 = new TLegend(0.6,0.6,0.89,0.89);
  le1->SetFillColor(10);
  le1->AddEntry(g10,"Thomas: 5   mm","p");
  le1->AddEntry(g20,"Thomas: 7.5 mm","p");
  le1->AddEntry(g30,"Thomas: 10  mm","p");
  le1->AddEntry(g40,"Thomas: 15  mm","p");
  le1->AddEntry(g50,"Thomas: 20  mm","p");
  le1->AddEntry(g60,"Thomas: 30  mm","p");
  // le1->AddEntry(g70,"Thomas:45  mm","p");
  le1->AddEntry(g80,"Thomas: 60  mm","p");
  le1->AddEntry(g90,"Yichen: 6.5 mm","p");

  le1->AddEntry(g100,"Cal. 89 K","l");
  le1->AddEntry(g200,"Cal. 86 K","l");
  le1->AddEntry(g300,"Cal. 92 K","l");
  le1->Draw();


 
  
  // c2->cd(2);
  
  // TH2F *h2 = new TH2F("h2","h2",100,0.09,4.0,100,0.00,2.);
  // h2->Draw();
  // h2->SetTitle("Applied Correction");
  // h2->SetXTitle("E field (kV/cm)");
  // h2->SetYTitle("#mu (cm^{2}/V/s)");
  // h2->GetYaxis()->SetTitleOffset(1.3);
  // TGraph *g10 = new TGraph(ncount[0],mE[0],norm_drift_time[0]);
  // TGraph *g20 = new TGraph(ncount[1],mE[1],norm_drift_time[1]);
  // TGraph *g30 = new TGraph(ncount[2],mE[2],norm_drift_time[2]);
  // TGraph *g40 = new TGraph(ncount[3],mE[3],norm_drift_time[3]);
  // TGraph *g50 = new TGraph(ncount[4],mE[4],norm_drift_time[4]);
  // TGraph *g60 = new TGraph(ncount[5],mE[5],norm_drift_time[5]);
  // TGraph *g70 = new TGraph(ncount[6],mE[6],norm_drift_time[6]);
  // TGraph *g80 = new TGraph(ncount[7],mE[7],norm_drift_time[7]);
  // TGraph *g90 = new TGraph(ncount[8],mE[8],norm_drift_time[8]);

  // g10->Draw("*same");
  // g10->SetMarkerStyle(20);
  // g10->SetMarkerColor(1);
  // g10->SetMarkerSize(1.0);
  // g20->Draw("*same");
  // g20->SetMarkerStyle(21);
  // g20->SetMarkerColor(3);
  // g30->Draw("*same");
  // g30->SetMarkerStyle(22);
  // g30->SetMarkerColor(4);
  // g40->Draw("*same");
  // g40->SetMarkerStyle(23);
  // g40->SetMarkerColor(6);
  
  // g50->Draw("*same");
  // g50->SetMarkerStyle(24);
  // g50->SetMarkerColor(1);

  // g60->Draw("*same");
  // g60->SetMarkerStyle(24);
  // g60->SetMarkerColor(8);
  // // g70->Draw("*same");
  // // g70->SetMarkerStyle(25);
  // // g70->SetMarkerColor(2);
  
  // g80->Draw("*same");
  // g80->SetMarkerStyle(25);
  // g80->SetMarkerColor(2);

  // g90->Draw("*same");
  // g90->SetMarkerStyle(25);
  // g90->SetMarkerColor(4);

}
