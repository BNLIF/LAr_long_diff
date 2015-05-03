#include "LAr.C"

void plot_result_sar(){
  Double_t cons[2],cons_err[2];
  Double_t slope[9],slope_err[9];
  ifstream infile1("bestfit_sar.dat");
  for (Int_t i=0;i!=2;i++){
    infile1 >> cons[i] >> cons_err[i];
    cons[i]  =fabs(cons[i]);
  }

 
  

  Double_t dis[7]={5,10,15,20,30,45,60};
  Double_t E[13],E1[13];
  
  ifstream infile("fit_result_sar.dat");
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

  for (Int_t i=0;i!=8;i++){
    infile1 >> slope[i] >> slope_err[i];
    slope[i] = slope[i];
    slope_err[i] = slope_err[i];
  }
  
  
  

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  
  c1->SetLogx(1);
  c1->SetLogy(1);
  TGraphErrors *g2 = new TGraphErrors(8,E,slope,0,slope_err);
  g2->Draw("A*");
  
  g2->GetXaxis()->SetRangeUser(0.07,3.5);
  g2->GetYaxis()->SetRangeUser(0.1,5.0);
  
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(20);
  g2->SetTitle("Electron Energy #epsilon_{L} (eV)");
  g2->GetXaxis()->SetTitle("E field (kV/cm)");
  g2->SetLineWidth(2.5);
 
  
  TLegend *le1 = new TLegend(0.15,0.6,0.5,0.89);
  le1->SetFillColor(10);
  le1->AddEntry(g2,"SAr (Long)","pl");
  
  le1->Draw();

  TFile *file = new TFile("sar.root","RECREATE");
  g2->Write("sar");
  file->Write();
  file->Close();
  
}
