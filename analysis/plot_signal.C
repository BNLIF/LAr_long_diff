void plot_signal(){
  TFile *file = new TFile("thomas.root");
  TTree *T = (TTree*)file->Get("T");
  Int_t fileno;
  Double_t E,dis;
  T->SetBranchAddress("fileno",&fileno);
  T->SetBranchAddress("E",&E);
  T->SetBranchAddress("dis",&dis);
  
  Int_t fno[100];
  Double_t Eno[100];
  Int_t no = 0;
  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    if (i<100&&dis==10){
      fno[no] = fileno;
      Eno[no] = E;
      no++;
    }
  }
  //cout << no << endl;
    
  TH1D **hh = new TH1D*[20];
  TF1 *f1 = new TF1("f1","pol0");
  Double_t par[10];
  for (Int_t i=0;i!=no;i++){
    //    cout << fno[i] << endl;
    hh[i] = (TH1D*)file->Get(Form("S1_%d",fno[i]));
    hh[i]->Fit(f1,"Q0","",-3500,-500);
    f1->GetParameters(par);
    cout << par[0] << endl;
    for (Int_t j=0;j!=hh[i]->GetNbinsX();j++){
      Double_t content = hh[i]->GetBinContent(j+1);
      content -= par[0];
      hh[i]->SetBinContent(j+1,content);
    }

    if (i==1){
      hh[i]->Draw();
    }else if (i>0){
      hh[i]->Draw("same");
    }
    hh[i]->SetLineColor(i+1);
    hh[i]->SetLineWidth(2.5);
  }
  hh[1]->GetYaxis()->SetRangeUser(-0.003,0.03);
  hh[1]->SetTitle("20 mm Drift Distance");
  hh[1]->SetXTitle("Time (ns)");
  hh[1]->SetYTitle("Signal (mV)");
  hh[1]->GetXaxis()->SetNdivisions(506);
  hh[1]->GetYaxis()->SetNdivisions(506);
  hh[1]->GetXaxis()->SetTitleSize(0.07);
  hh[1]->GetXaxis()->SetLabelSize(0.07);
  hh[1]->GetYaxis()->SetTitleSize(0.07);
  hh[1]->GetYaxis()->SetLabelSize(0.07);
  


  hh[9]->SetLineColor(1);
  TLegend *le1 = new TLegend(0.6,0.15,0.89,0.89);
  le1->SetFillColor(10);
  for (Int_t i=1;i!=no;i++){
    le1->AddEntry(hh[i],Form("E = %3.2f kV/cm",Eno[i]),"l");
  }
  le1->Draw();
}
