#include "TH1.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"
#include "TStopwatch.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TRotation.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "TCut.h"
#include "TRandom.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include "TGraph.h"

using std::vector;
using namespace std;

class MyFCN : public ROOT::Minuit2::FCNBase { 
  
public: 
  double Up() const { return 1.; }
  
  MyFCN() {
  }

  ~MyFCN(){
  }


  void init(){
    ifstream infile("fit_result.dat");
    npoint = 0;
    while(!infile.eof()){
      infile >> num[npoint] >> E[npoint]
	     >> dis[npoint] >> time[npoint] >> time_err[npoint]
	     >> dtime[npoint] >> dtime_err[npoint];
      dtime_err[npoint] *= dtime[npoint];
      npoint ++;
    }
    npoint --;
    infile.close();
  }

  double operator() (const std::vector<double> & x) const {
    Double_t results = 0;
    
    Double_t cons[7];
    Double_t slope[40];
    
    for (Int_t i=0;i!=7;i++){
      cons[i] = x[i];
    }
    for (Int_t i=7;i!=7+40;i++){
      slope[i-7] = x[i]/100.; 
    }

    for (Int_t i=0;i!=npoint;i++){
      Double_t con=0;
      if (dis[i]==5){
	con = cons[0];
      }else if (dis[i]==7.5){
	con = cons[1];
      }else if (dis[i]==10.){
	con = cons[2];
      }else if (dis[i]==15.){
	con = cons[3];
      }else if (dis[i]==20.){
	con = cons[4];
      }else if (dis[i]==30.){
	con = cons[5];
      }else if (dis[i]==60.){
	con = cons[6];
      }

      Double_t diff = pow(time[i],2) -(con*con + 2*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i]);
      
      Double_t diff_err2 = pow(2*time[i],2) * (time_err[i]*time_err[i] + pow(time[i]*0.09,2))
	+ pow(2*slope[num[i]]/dis[i]/E[i],2)*pow(2*dtime[i]*dtime_err[i],2);
      results += diff * diff / diff_err2;
      
    }


    return results;
  }

   double getchi2(double *x) {
    Double_t results = 0;
    
    Double_t cons[7];
    Double_t slope[40];
    
    for (Int_t i=0;i!=7;i++){
      cons[i] = x[i];
    }
    for (Int_t i=7;i!=7+40;i++){
      slope[i-7] = x[i]/100.; 
    }

    for (Int_t i=0;i!=npoint;i++){
      Double_t con=0;
      if (dis[i]==5){
	con = cons[0];
      }else if (dis[i]==7.5){
	con = cons[1];
      }else if (dis[i]==10.){
	con = cons[2];
      }else if (dis[i]==15.){
	con = cons[3];
      }else if (dis[i]==20.){
	con = cons[4];
      }else if (dis[i]==30.){
	con = cons[5];
      }else if (dis[i]==60.){
	con = cons[6];
      }
      
      //cout << num[i] << " " << dis[i] << " " << E[i] << " " << con << " " << slope[num[i]] << endl;
      
      Double_t diff = pow(time[i],2) -(con*con + 2*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i]);
      
      Double_t diff_err2 = pow(2*time[i],2) * (time_err[i]*time_err[i] + pow(time[i]*0.09,2))
	+ pow(2*slope[num[i]]/dis[i]/E[i],2)*pow(2*dtime[i]*dtime_err[i],2);
      
      results += diff * diff / diff_err2;
      
    }
    
    return results;
   }

private:
  int npoint;
  int num[100];
  Double_t E[100];
  Double_t dis[100];
  Double_t time[100],time_err[100];
  Double_t dtime[100],dtime_err[100];
};




int testMinimize() { 
  

  TStopwatch w; 
  
  double par[300];
  for (Int_t i=1;i!=300;i++){
    par[i] = 0.;
  }
  // 7 parameters for the constant background
  // 15 parameters for 15 E-field settings

  TFitterMinuit *minuit_shape = new TFitterMinuit(); 
  MyFCN fcn_shape;
  fcn_shape.init();
  minuit_shape->SetMinuitFCN(&fcn_shape);

  
  par[0] = 50;
  par[1] = 70;
  par[2] = 60;
  par[3] = 40;
  par[4] = 80;
  par[5] = 50;
  par[6] = 50;

  // par[7] = 124577;
  // par[8] = 21884.1;
  // par[9] = 15511.6;
  // par[10] = 8425.29;
  // par[11] = 2321.15;
  // par[12] = 1678.33;

  
  // par[13] = 0.0150277;
  // par[14] = 0.0104349;
  // par[15] = 0.0155893;
  // par[16] = 0.0197897;
  // par[17] = 0.0224081;
  // par[18] = 0.0338801;
  // par[19] = 0.0281255;
  // par[20] = 0.0336592;
  // par[21] = 0.0426773;
  // par[22] = 0.0366487;
  // par[23] = 0.0479296;
  // par[24] = 0.0364158;
  // par[25] = 1.40037e-05;
  
  
  for (Int_t i=0;i!=7;i++){
    minuit_shape->SetParameter(i,Form("con%d",i),par[i],0.1,0,150.);
  }
  for (Int_t i=7;i!=7+40;i++){
    minuit_shape->SetParameter(i,Form("E%d",i),par[i],0.0001,0,0.2);
  }
  for (Int_t i=7;i!=7+13;i++){
    minuit_shape->FixParameter(i);
  }
  
  minuit_shape->SetPrintLevel(3);
  minuit_shape->CreateMinimizer(TFitterMinuit::kCombined);
  w.Start();
  int iret = minuit_shape->Minimize(0,1.0);
  w.Stop();
  
  ofstream outfile("bestfit.dat");
  for (Int_t i=0;i!=7+40;i++){
    par[i] = minuit_shape->GetParameter(i);
    outfile << par[i] << " " << minuit_shape->GetParError(i) << endl;
  } 
  fcn_shape.getchi2(par);

  
  
  
  std::cout << "\nTime: \t" << w.RealTime() << " , " << w.CpuTime() << std::endl; 

  return 0;
}



int main(Int_t argc, char *argv[]){
  if (argc==1){
    cout << "Command Line: ./fit -a[x_init] -b[neve]" << endl;
    return 0;
  }else{
    
    
    int iret = testMinimize();
    
    if (iret != 0) { 
      std::cerr << "ERROR: Minimize test failed !" << std::endl;
      return iret;
    }
  }
}
