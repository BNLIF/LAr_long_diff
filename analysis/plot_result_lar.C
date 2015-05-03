#include "LAr.C"

void plot_result_lar(){
  Double_t cons[8],cons_err[8];
  Double_t slope[29],slope_err[29];
  ifstream infile1("bestfit_lar.dat");
  for (Int_t i=0;i!=8;i++){
    infile1 >> cons[i] >> cons_err[i];
    cons[i]  =fabs(cons[i]);
  }

 
  

  Double_t dis[8]={5,7.5,10,15,20,30,45,60};
  Double_t E[29],E1[29];
  
  ifstream infile("fit_result_lar.dat");
  Int_t num;
  Double_t EE;
  Double_t temp;
  while(!infile.eof()){
    infile >> num >> EE 
	   >> temp >> temp >> temp
	   >>temp >> temp;
    E[num] = EE;
    E1[num] = EE;
   
  }

  LAr argon;
  Double_t temperature = 87;
  Double_t calc[29];
  Double_t calcT[29];

  for (Int_t i=0;i!=29;i++){
    infile1 >> slope[i] >> slope_err[i];

    Double_t vel = argon.vDrift(temperature,E[i]);
    calc[i] = argon.Diffusion(temperature,E[i],2);
    calcT[i] = argon.Diffusion(temperature,E[i],1);
    slope[i] = slope[i];
    slope_err[i] = slope_err[i];
    
    //cout << E[i] << " " << slope[i] << endl;
  }
  slope[25]=0;
  Double_t a1,a2,a3;
  a1 = 1./slope_err[11];
  a2 = 1./slope_err[12];
  E[12] = (a1*a1*E[11]+a2*a2*E[12])/(a1*a1+a2*a2);
  slope[12] = (a1*a1*slope[11]+a2*a2*slope[12])/(a1*a1+a2*a2);
  slope_err[12] = 1./sqrt(pow(a1,2)+pow(a2,2));
  slope[11] = 0;

  a1 = 1./slope_err[13];
  a2 = 1./slope_err[14];
  E[14] = (a1*a1*E[13]+a2*a2*E[14])/(a1*a1+a2*a2);
  slope[14] = (a1*a1*slope[13]+a2*a2*slope[14])/(a1*a1+a2*a2);
  slope_err[14] = 1./sqrt(pow(a1,2)+pow(a2,2));
  slope[13] = 0;

  a1 = 1./slope_err[15];
  a2 = 1./slope_err[16];
  E[16] = (a1*a1*E[15]+a2*a2*E[16])/(a1*a1+a2*a2);
  slope[16] = (a1*a1*slope[15]+a2*a2*slope[16])/(a1*a1+a2*a2);
  slope_err[16] = 1./sqrt(pow(a1,2)+pow(a2,2));
  slope[15] = 0;

  a1 = 1./slope_err[17];
  a2 = 1./slope_err[18];
  E[18] = (a1*a1*E[17]+a2*a2*E[18])/(a1*a1+a2*a2);
  slope[18] = (a1*a1*slope[17]+a2*a2*slope[18])/(a1*a1+a2*a2);
  slope_err[18] = 1./sqrt(pow(a1,2)+pow(a2,2));
  slope[17] = 0;


  a1 = 1./slope_err[19];
  a2 = 1./slope_err[20];
  a3 = 1./slope_err[21];
  E[21] = (a1*a1*E[19]+a2*a2*E[20] + a3*a3*E[21])/(a1*a1+a2*a2+a3*a3);
  slope[21] = (a1*a1*slope[19]+a2*a2*slope[20]+a3*a3*slope[21])/(a1*a1+a2*a2+a3*a3);
  slope_err[21] = 1./sqrt(pow(a1,2)+pow(a2,2)+pow(a3,2));
  slope[19] = 0;
  slope[20] = 0;

  a1 = 1./slope_err[22];
  a2 = 1./slope_err[23];
  E[23] = (a1*a1*E[22]+a2*a2*E[23])/(a1*a1+a2*a2);
  slope[23] = (a1*a1*slope[22]+a2*a2*slope[23])/(a1*a1+a2*a2);
  slope_err[23] = 1./sqrt(pow(a1,2)+pow(a2,2));
  slope[22] = 0;
 

  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->Divide(3,1);
  // c1->cd(1);
  // TGraphErrors *g1 = new TGraphErrors(8,dis,cons,0,cons_err);
  // g1->Draw("A*");
  // g1->GetXaxis()->SetTitle("Drift Distance (mm)");
  // g1->GetYaxis()->SetTitle("c (ns)");
  // g1->SetTitle();
  // g1->SetMarkerColor(2);
  // g1->SetMarkerStyle(20);
  // g1->SetLineWidth(2.5);
  // g1->GetYaxis()->SetRangeUser(0.,120.);
  // g1->GetYaxis()->SetTitleOffset(1.2);
 


  c1->cd(1);
  
  TFile *file1 = new TFile("gar.root");
  TGraphErrors *g1 = (TGraphErrors*)file1->Get("gar");
  g1->Draw("A*");
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetTitle("GAr");
  g1->GetYaxis()->SetTitle("Electron Energy #epsilon_{L} eV");
  g1->GetXaxis()->SetTitle("E/P (V/cm/torr)");

  // c1->cd(2);
  c1->cd(2);
  c1_2->SetLogx(1);
  //  c1_2->SetLogy(1);
  TGraphErrors *g2 = new TGraphErrors(29,E,slope,0,slope_err);
  g2->Draw("A*");
  TGraph *g3 = new TGraph(29,E1,calc);
  g3->Draw("Lsame");
  g2->GetXaxis()->SetRangeUser(0.07,3.5);
  g2->GetYaxis()->SetRangeUser(0.005,0.1);
  TGraph *g4 = new TGraph(29,E1,calcT);
  g4->Draw("Lsame");
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(20);
  g2->SetTitle("LAr");
  g2->GetXaxis()->SetTitle("E field (kV/cm)");
  g2->SetLineWidth(2.5);
  g3->SetLineWidth(2.5);
  g3->SetLineColor(6);
  g4->SetLineWidth(2.5);
  g4->SetLineColor(4);
  
  TLegend *le1 = new TLegend(0.15,0.6,0.5,0.89);
  le1->SetFillColor(10);
  le1->AddEntry(g2,"This work (Long)","pl");
  le1->AddEntry(g3,"Long. Cal.","pl");
  le1->AddEntry(g4,"Tran. Cal.","pl");
  le1->Draw();

  c1->cd(3);
  TFile *file1 = new TFile("sar.root");
  TGraphErrors *g1 = (TGraphErrors*)file1->Get("sar");
  g1->Draw("A*");
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetTitle("SAr");
  g1->GetYaxis()->SetTitle("");
  g1->GetXaxis()->SetTitle("E (kV/cm)");
  
}
