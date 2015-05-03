#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

Double_t riseFunc(Double_t *x, Double_t*par){
    Double_t a  = par[0];
    Double_t tI = par[1];
    Double_t tD = par[2];
    Double_t c  = par[3];
    Double_t s  = par[4];
    Double_t bl = par[5];
    Double_t flag = par[6];
 
    Double_t pos1 = (s*s-2*(x[0]-c)*tD)/(2*tD*tD);
    Double_t pos2 = (s*s-2*(x[0]-c)*tI)/(2*tI*tI);
    Double_t pos3 = (-s+(x[0]-c)*tD/s)/TMath::Sqrt(2)/tD;
    Double_t pos4 = (-s+(x[0]-c)*tI/s)/TMath::Sqrt(2)/tI;

    double f =
      TMath::Exp(pos1)*(1+TMath::Erf(pos3))
      + TMath::Exp(pos2)*(-2+TMath::Erfc(pos4))
      ;

    TF1 ff("test",Form("TMath::Exp(-x*x)*x/sqrt(x*x+%f)",pos2),8500,9500);
    ROOT::Math::WrappedTF1 wf1(ff);
    ROOT::Math::GaussIntegrator ig;
    ig.SetFunction(wf1);
    ig.SetRelTolerance(0.001);
    Double_t f1;
    f1 = -2./sqrt(3.1415926)*ig.Integral(sqrt(pos4*pos4-pos2),100);
    if (pos4>0){
      f1 += 2*(-2./sqrt(3.1415926))*ig.Integral(sqrt(-pos2),sqrt(pos4*pos4-pos2));
    }
    f1 += TMath::Exp(pos1)*(1+TMath::Erf(pos3));

    if (flag==1) { return pos1; }
    else if (flag==2) { return pos2;}
    else if (flag==3) { return pos3;}
    else if (flag==4) { return pos4;}
    else if (flag==5) { return TMath::Exp(pos1);}
    else if (flag==6) { return TMath::Exp(pos2);}
    else if (flag==7) { return (1+TMath::Erf(pos3)) ;}
    else if (flag==8) { return fabs((-2+TMath::Erfc(pos4)));}
    else if (flag==9) { return f;}
    else if (flag==10) { return f1;}

    
}



void test(){
  Double_t a = 1;  // normalization
  Double_t tI = 2.29765e+01;  //inpulse function time
  Double_t tD = 1.10072e5; // decay function time
  Double_t c  = 9.00296e3; // where the time is ... 
  Double_t s  = 1.92314e2; //sigma to be fitted
  Double_t bl = 0;  // baseline 
  
  Double_t par[10];
  par[0] = a;
  par[1] = tI;
  par[2] = tD;
  par[3] = c;
  par[4] = s;
  par[5] = bl;
  par[6] = 1;
  
  TF1 *f1 = new TF1("f1",riseFunc,8500,9500,7);
  TF1 *f2 = new TF1("f2",riseFunc,8500,9500,7);
  TF1 *f3 = new TF1("f3",riseFunc,8500,9500,7);
  TF1 *f4 = new TF1("f4",riseFunc,8500,9500,7);
  TF1 *f5 = new TF1("f5",riseFunc,8500,9500,7);
  TF1 *f6 = new TF1("f6",riseFunc,8500,9500,7);
  TF1 *f7 = new TF1("f7",riseFunc,8500,9500,7);
  TF1 *f8 = new TF1("f8",riseFunc,8500,9500,7);
  TF1 *f9 = new TF1("f9",riseFunc,8500,9500,7);
  TF1 *f10 = new TF1("f10",riseFunc,8500,9500,7);
  f1->SetParameters(par);
  par[6] = 2; f2->SetParameters(par);
  par[6] = 3; f3->SetParameters(par);
  par[6] = 4; f4->SetParameters(par);
  par[6] = 5; f5->SetParameters(par);
  par[6] = 6; f6->SetParameters(par);
  par[6] = 7; f7->SetParameters(par);
  par[6] = 8; f8->SetParameters(par);
  par[6] = 9; f9->SetParameters(par);
  par[6] = 10; f10->SetParameters(par);

  f10->Eval(8700) ;
  f10->Eval(9300) ;
  // cout << f5->Eval(9000) << endl;

  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(3,3);
  c1->cd(1);
  f1->Draw();
  c1->cd(2);
  f2->Draw();
  c1->cd(3);
  f3->Draw();

  c1->cd(4);
  f4->Draw();
  c1->cd(5);
  f5->Draw();
  c1->cd(6);
  f6->Draw();

  c1->cd(7);
  f7->Draw();
  c1->cd(8);
  f8->Draw();
  c1->cd(9);
  f9->Draw();
  f10->Draw("same");
  f10->SetLineColor(1);
}
