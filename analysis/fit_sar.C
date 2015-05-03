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
    ifstream infile("fit_result_sar.dat");
    npoint = 0;
    while(!infile.eof()){
      infile >> num[npoint] >> E[npoint]
	     >> dis[npoint] 
	     >> dtime[npoint] >> dtime_err[npoint]
	     >> time[npoint] >> time_err[npoint];
      
      // dtime[npoint] = sqrt(dtime[npoint]*dtime[npoint]-32*32-20*20);
      // dtime_err[npoint] = sqrt(pow(dtime[npoint]*dtime_err[npoint],2)+pow(32*2,2)+pow(20*20,2))/dtime[npoint];
      
      npoint ++;
    }
    npoint --;
    cout << npoint << endl;
    infile.close();
  }

  double operator() (const std::vector<double> & x) const {
    Double_t results = 0;
    
    Double_t cons[2];
    Double_t slope[9];
    
    for (Int_t i=0;i!=2;i++){
      cons[i] = x[i];
    }
    for (Int_t i=2;i!=2+8;i++){
      slope[i-2] = x[i]/100.; 
    }

    for (Int_t i=0;i!=npoint;i++){
      Double_t con=0;
      if (dis[i]==15){
      	con = cons[0];
      }else if (dis[i]==30.){
      	con = cons[1];
      }

      Double_t diff = pow(time[i],2) -(con*con + 2.*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i]);
      
      Double_t diff_err2 = pow(2*time[i],2) * (time_err[i]*time_err[i] + pow(time[i]*0.0,2))
	+ pow(2*slope[num[i]]/dis[i]/E[i],2)*pow(2*dtime[i]*dtime_err[i],2)
	+ pow(2.*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i],2)*pow(0.1/dis[i],2)
	+ pow(2.*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i],2)*pow(0.01,2);
      results += diff * diff / diff_err2 + pow(con - 38,2)/pow(11.,2);
      
    }


    return results;
  }

   double getchi2(double *x) {
    Double_t results = 0;
    Double_t cons[2];
    Double_t slope[9];
    
    for (Int_t i=0;i!=2;i++){
      cons[i] = x[i];
    }
    for (Int_t i=2;i!=2+8;i++){
      slope[i-2] = x[i]/100.; 
    }

    for (Int_t i=0;i!=npoint;i++){
      Double_t con=0;
      // if (dis[i]==5){
      //con = cons[0];
      // }else if (dis[i]==10.){
      // 	con = cons[1];
      // }

      Double_t diff = pow(time[i],2) -(con*con + 2.*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i]);
      
      Double_t diff_err2 = pow(2*time[i],2) * (time_err[i]*time_err[i] + pow(time[i]*0.0,2))
	+ pow(2*slope[num[i]]/dis[i]/E[i],2)*pow(2*dtime[i]*dtime_err[i],2)
	+ pow(2.*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i],2)*pow(0.1/dis[i],2)
	+ pow(2.*slope[num[i]]*pow(dtime[i],2)/dis[i]/E[i],2)*pow(0.01,2);
      results += diff * diff / diff_err2 + pow(con - 38,2)/pow(11.,2);
      
      cout << diff << " " << sqrt(diff_err2) << endl;
    }
    cout << results << endl;



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
  // 13 parameters for 13 E-field settings

  TFitterMinuit *minuit_shape = new TFitterMinuit(); 
  MyFCN fcn_shape;
  fcn_shape.init();
  minuit_shape->SetMinuitFCN(&fcn_shape);

  
  // par[0] = 50;
  // par[1] = 70;
  // par[2] = 60;
  // par[3] = 40;
  // par[4] = 80;
  // par[5] = 50;
  // par[6] = 50;
 
  
  
  
  for (Int_t i=0;i!=2;i++){
    minuit_shape->SetParameter(i,Form("con%d",i),par[i],0.1,-0,500.);
    //minuit_shape->FixParameter(i);
  }
  for (Int_t i=2;i!=2+8;i++){
    minuit_shape->SetParameter(i,Form("E%d",i),par[i],0.0001,0,5.0);
  }
  
  minuit_shape->SetPrintLevel(3);
  minuit_shape->CreateMinimizer(TFitterMinuit::kCombined);
  w.Start();
  int iret = minuit_shape->Minimize(0,1.0);
  w.Stop();
  
  ofstream outfile("bestfit_sar.dat");
  for (Int_t i=0;i!=2+8;i++){
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
