#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"

void test1(){
  Double_t pos2 = 10;
  Double_t pos4 = 3;

  TF1 ff("test",Form("TMath::Exp(-x*x)*x/sqrt(x*x+%f)",pos2),8500,9500);
  ROOT::Math::WrappedTF1 wf1(ff);
  ROOT::Math::GaussIntegrator ig;
  ig.SetFunction(wf1);
  ig.SetRelTolerance(0.001);
  Double_t f1 = ig.Integral(pos4*pos4-pos2,100);
  cout << f1 << endl;
}
