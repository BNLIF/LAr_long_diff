//fit parameters:
//                  par0 = rise edge amplitude;
//                  par1 = Impulse reponse time constant 28.32ns for our preamp;
//                  par2 = Decay time constant of the preamp;
//                  par3 = drift time;
//                  par4 = rise time;
//                  par5 = base line shift;

#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"
#include "TStopwatch.h"

Double_t riseFunc(Double_t *x, Double_t*par){
    Double_t a  = par[0];
    Double_t tI = par[1];
    Double_t tD = par[2];
    Double_t c  = par[3];
    Double_t s  = par[4];
    Double_t bl = par[5];
    Double_t min = par[6];
    Double_t max = par[7];
    
    Double_t pos1 = (s*s-2*(x[0]-c)*tD)/(2*tD*tD);
    Double_t pos2 = (s*s-2*(x[0]-c)*tI)/(2*tI*tI);
    Double_t pos3 = (-s+(x[0]-c)*tD/s)/TMath::Sqrt(2)/tD;
    Double_t pos4 = (-s+(x[0]-c)*tI/s)/TMath::Sqrt(2)/tI;

    
    // TF1 ff("test",Form("TMath::Exp(-x*x)*x/sqrt(x*x+%f)",pos2),8500,9500);
    // ROOT::Math::WrappedTF1 wf1(ff);
    // ROOT::Math::GaussIntegrator ig;
    // ig.SetFunction(wf1);
    // ig.SetRelTolerance(0.001);
    Double_t f1=0;
    // f1 = -2./sqrt(3.1415926)*ig.Integral(sqrt(pos4*pos4-pos2),100);
    // if (pos4>0){
    //   f1 += 2*(-2./sqrt(3.1415926))*ig.Integral(sqrt(-pos2),sqrt(pos4*pos4-pos2));
    // }
    f1 += TMath::Exp(pos1)*(1+TMath::Erf(pos3));

    double f = bl+a/2*f1;


    return f; 
}


void fit_test1(Int_t num = 0, Int_t flag = 0, Int_t fit_min = 200, Int_t fit_max = 6000, Double_t fit_time = 950, Double_t fit_diff = 50, Double_t height = 0.01){
  TFile *file = new TFile("thomas.root");
  TString name;
  name.Form("S1_%d",num);
  TH1D *h1 = (TH1D*)file->Get(name);
  h1->Draw();


  Double_t max = -100;
  Double_t min = 100;
  for (Int_t i=0;i!=h1->GetNbinsX()-1;i++){
    Double_t prev_content = h1->GetBinContent(i+1);
    Double_t content = h1->GetBinContent(i+100);
    Double_t diff = content - prev_content;
    if (diff < min) min = diff;
    if (diff > max) max = diff;
  }
  TH1F *h2 = new TH1F("h2","h2",100,min,-min);
  
  for (Int_t i=0;i!=h1->GetNbinsX()-1;i++){
    Double_t prev_content = h1->GetBinContent(i+1);
    Double_t content = h1->GetBinContent(i+100);
    Double_t diff = content - prev_content;
    h2->Fill(diff);
  }
  TH1F *h3 = new TH1F("h3","h3",50,-h2->GetRMS()*3,h2->GetRMS()*3);
 
  for (Int_t i=0;i!=h1->GetNbinsX()-1;i++){
    Double_t prev_content = h1->GetBinContent(i+1);
    Double_t content = h1->GetBinContent(i+100);
    Double_t diff = content - prev_content;
    if (fabs(diff)<h2->GetRMS()*3)
      h3->Fill(diff);
  }
  h3->Draw();
  // h3->Fit("gaus");
  Double_t err = h3->GetRMS();

  cout << h2->GetRMS() << endl;
  cout << err << endl;


  // when apply error, use a sqrt(2.)
  TH1D *h4 = (TH1D*)h1->Clone("h4");
  for (Int_t i=0;i!=h4->GetNbinsX();i++){
    h4->SetBinError(i+1,err/sqrt(2.));
  }
  h4->Draw();

  
  TFile *file = new TFile("result_thomas_new.root");
  TTree *T = (TTree*)file->Get("T");
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

  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (fileno==num) break;
  }

  cout << "Drift Time: " << m_drift_time <<  endl;
  //cout << max << endl;

  TLine *l1 = new TLine(m_drift_time,-1,m_drift_time,1);
  l1->Draw();
  l1->SetLineColor(2);

  if (flag==1){
    TF1 *f1 = new TF1("f1",riseFunc,h4->GetBinCenter(1),h4->GetBinCenter(h4->GetNbinsX()),6);
    Double_t par[6]={0.001,  //height
		     22.9765, //fixed
		     250000, // decay time
		     60, // drift time
		     10, //results
		     0}; //baseline
    
    //par[0] = 1;

    // par[0] = max;
    // par[3] = fabs(m_drift_time);
    // par[4] = fabs(m_diff_time);

    par[0] = height;
    par[3] = fit_time;
    par[4] = fit_diff;

    // par[0] = 0.2e-3;
    // par[3] = 260000;
    // par[4] = 2400;
    par[5] = -0.2e-3;
    
    f1->SetParameters(par);
    f1->FixParameter(1,par[1]);
    //f1->Draw("same");
    
    
    
    TStopwatch w; 
    w.Start();
    h4->Fit("f1","R","",fit_min,fit_max);
    w.Stop();
    //std::cout << "\nTime: \t" << w.RealTime() << " , " << w.CpuTime() << std::endl; 
    cout << num << " " << fit_min << " " << fit_max << " " << f1->GetParameter(4) << " " << f1->GetParError(4) << " " << f1->GetChisquare()/f1->GetNDF() << " " 
	 << f1->GetParameter(3) << " " << f1->GetParError(3)	 << endl;
  }else{
    TLine *l1 = new TLine(fit_min,-1,fit_min,1);
    TLine *l2 = new TLine(fit_max,-1,fit_max,1);
    l1->Draw("same");
    l2->Draw("same");
  }

  

}

