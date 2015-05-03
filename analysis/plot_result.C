#include "LAr.C"

void plot_result(){
  Double_t cons[7],cons_err[7];
  Double_t slope[40],slope_err[40];
  ifstream infile1("bestfit.dat");
  for (Int_t i=0;i!=7;i++){
    infile1 >> cons[i] >> cons_err[i];
  }
  Double_t dis[7]={5,7.5,10,15,20,30,60};
  Double_t E[40]={0.01,0.01,0.01,0.01,0.01,
		  0.01,0.01,0.01,0.01,0.01,
		  0.01,0.01,0.01,0.1,0.13296,
		  0.150231,0.200244,0.266907,0.299419,0.399098,
		  0.450941,0.498655,0.531961,0.601063,0.664661,
		  0.801161,0.931623,1.00101,1.09901,1.19795,
		  1.40307,1.46487,1.49679,1.99508,2.12834,
		  2.49276,2.7964,3.00469,3.46897,4.00498
		  };

  LAr argon;
  Double_t temperature = 89;
  Double_t calc[40];
  Double_t calcT[40];

  for (Int_t i=0;i!=40;i++){
    infile1 >> slope[i] >> slope_err[i];

    Double_t vel = argon.vDrift(temperature,E[i]);
    calc[i] = argon.Diffusion(temperature,E[i],2);
    calcT[i] = argon.Diffusion(temperature,E[i],1);
    slope[i] = slope[i];
    slope_err[i] = slope_err[i];
    
  }
  
 

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  TGraphErrors *g1 = new TGraphErrors(7,dis,cons,0,cons_err);
  g1->Draw("A*");
  g1->GetXaxis()->SetTitle("Drift Distance (mm)");
  g1->GetYaxis()->SetTitle("#sigma_{const} (ns)");
  g1->SetTitle();
  g1->SetMarkerColor(2);
  g1->SetMarkerStyle(20);
  g1->SetLineWidth(2.5);
 
  
  c1->cd(2);
  c1_2->SetLogx(1);
  c1_2->SetLogy(1);
  TGraphErrors *g2 = new TGraphErrors(40,E,slope,0,slope_err);
  g2->Draw("A*");
  TGraph *g3 = new TGraph(40,E,calc);
  g3->Draw("Lsame");
  g2->GetXaxis()->SetRangeUser(0.09,3.5);
  g2->GetYaxis()->SetRangeUser(0.007,0.1);
  TGraph *g4 = new TGraph(40,E,calcT);
  g4->Draw("Lsame");
  g2->SetMarkerColor(2);
  g2->SetMarkerStyle(20);
  g2->SetTitle("Electron Energy (eV)");
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
  
}
