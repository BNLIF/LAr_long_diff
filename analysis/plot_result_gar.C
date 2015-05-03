#include "LAr.C"

void plot_result_gar(){
  Double_t cons[7],cons_err[7];
  Double_t slope[13],slope_err[13];
  ifstream infile1("bestfit_gar.dat");
  for (Int_t i=0;i!=7;i++){
    infile1 >> cons[i] >> cons_err[i];
    cons[i]  =fabs(cons[i]);
  }

 
  

  Double_t dis[7]={5,10,15,20,30,45,60};
  Double_t E[13],E1[13];
  
  ifstream infile("fit_result_gar.dat");
  Int_t num;
  Double_t EE;
  Double_t temp;
  while(!infile.eof()){
    infile >> num >> EE 
	   >> temp >> temp >> temp
	   >>temp >> temp;
    E[num] = EE*1000./1134.8;
    E1[num] = EE*1000./1134.8;
   
  }

  LAr argon;
  Double_t temperature = 87;
  Double_t calc[29];
  Double_t calcT[29];

  for (Int_t i=0;i!=13;i++){
    infile1 >> slope[i] >> slope_err[i];


    slope[i] = slope[i];
    slope_err[i] = slope_err[i];
    
  }
  
 

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  
  c1->SetLogx(1);
  c1->SetLogy(1);
  TGraphErrors *g2 = new TGraphErrors(13,E,slope,0,slope_err);
  g2->Draw("A*");
  
  g2->GetXaxis()->SetRangeUser(0.07,3.5);
  g2->GetYaxis()->SetRangeUser(0.1,5.0);
  
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(20);
  g2->SetTitle("Electron Energy #epsilon_{L} (eV)");
  g2->GetXaxis()->SetTitle("E field (kV/cm)");
  g2->SetLineWidth(2.5);
 
  TFile *file = new TFile("gar.root","RECREATE");
  g2->Write("gar");
  file->Write();
  file->Close();
  // TGraphErrors *g3 = (TGraphErrors*)file->Get("sar");

  // g3->Draw("*same");
  // g3->SetMarkerColor(4);
  // g3->SetMarkerStyle(21);
 
  // TLegend *le1 = new TLegend(0.6,0.15,0.89,0.5);
  // le1->SetFillColor(10);
  // le1->AddEntry(g2,"GAr (Long)","pl");
  // le1->AddEntry(g3,"SAr (Long)","pl");
  
  // le1->Draw();
  
}
