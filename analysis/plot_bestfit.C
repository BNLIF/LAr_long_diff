void plot_bestfit(){
  int npoint;
  int num[100];
  Double_t E[100];
  Double_t dis[100];
  Double_t time[100],time_err[100];
  Double_t dtime[100],dtime_err[100];

  ifstream infile("fit_result.dat");
  npoint = 0;
  while(!infile.eof()){
    infile >> num[npoint] >> E[npoint]
	   >> dis[npoint] >> time[npoint] >> time_err[npoint]
	   >> dtime[npoint] >> dtime_err[npoint];
    dtime_err[npoint] *= dtime[npoint];
    time_err[npoint] = sqrt(pow(time_err[npoint],2) + pow(time[npoint]*0.09,2));
    npoint ++;
  }
  npoint --;
  infile.close();

  Double_t cons[7],cons_err[7];
  Double_t slope[40],slope_err[40];
  ifstream infile1("bestfit.dat");
  for (Int_t i=0;i!=7;i++){
    infile1 >> cons[i] >> cons_err[i];
  }
  for (Int_t i=0;i!=40;i++){
    infile1 >> slope[i] >> slope_err[i];
  }
  
  TGraphErrors **gdata = new TGraphErrors*[7];
  TGraph **gfit = new TGraph*[7];
  Int_t ngpoint[7];
  for (Int_t i=0;i!=7;i++){
    gdata[i] = new TGraphErrors();
    gfit[i] = new TGraph();
    ngpoint[i] = 0;
  }
  
  for (Int_t i=0;i!=npoint;i++){
    Double_t con=0;
    Int_t ng = 0;
    if (dis[i]==5){
      con = cons[0];  ng = 0 ;
    }else if (dis[i]==7.5){
      con = cons[1];  ng = 1 ;
    }else if (dis[i]==10.){
      con = cons[2];  ng = 2 ;
    }else if (dis[i]==15.){
      con = cons[3];  ng = 3 ;
    }else if (dis[i]==20.){
      con = cons[4];  ng = 4 ;
    }else if (dis[i]==30.){
      con = cons[5];  ng = 5 ;
    }else if (dis[i]==60.){
      con = cons[6];  ng = 6 ;
    }
    
    gdata[ng]->SetPoint(ngpoint[ng],E[i],time[i]);
    gdata[ng]->SetPointError(ngpoint[ng],0,time_err[i]);
    gfit[ng]->SetPoint(ngpoint[ng],E[i],sqrt(con*con+2*slope[num[i]]/100.*pow(dtime[i],2)/dis[i]/E[i]));
    ngpoint[ng]++;

  }

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(3,3);
  c1->cd(1);
  
  Int_t i=0;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("5 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);

  c1->cd(2);
  i=1;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("7.5 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);

  c1->cd(3);
  i=2;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("10 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);

  c1->cd(4);
  i=3;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("15 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);

  c1->cd(5);
  i=4;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("20 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);

  c1->cd(6);
  i=5;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("30 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);
  
  c1->cd(7);
  i=6;
  gdata[i]->Draw("A*");
  gfit[i]->Draw("Lsame");
  gdata[i]->GetXaxis()->SetTitle("E field (kV/cm)");
  gdata[i]->GetYaxis()->SetTitle("T_{diff} (ns)");
  gdata[i]->SetTitle("60 mm Drift Distance");
  gdata[i]->SetMarkerStyle(20);
  gdata[i]->SetMarkerColor(2);
  gdata[i]->SetLineColor(1);
  gdata[i]->SetLineWidth(1.5);
  gfit[i]->SetLineColor(4);
  gfit[i]->SetLineWidth(1.5);
  gdata[i]->GetXaxis()->SetRangeUser(0.08,4);
}
