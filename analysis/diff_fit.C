#include<cmath>
//fit rising edge with our function

//fit parameters:
//                  par0 = rise edge amplitude;
//                  par1 = Impulse reponse time constant 28.32ns for our preamp;
//                  par2 = Decay time constant of the preamp;
//                  par3 = drift time;
//                  par4 = rise time;
//                  par5 = base line shift;
double riseFunc1(double *x, double*par){
    double a  = par[0];
    double tI = par[1];
    double tD = par[2];
    double c  = par[3];
    double s  = par[4];
    double bl = par[5];
    double f = bl+a/(2*(tD-tI))*tD*tI*((TMath::Exp((s*s-2*(x[0]-c)*tD)/(2*tD*tD))*(1+TMath::Erf((-s+(x[0]-c)*tD/s)/TMath::Sqrt(2)/tD)))+ (TMath::Exp((s*s-2*(x[0]-c)*tI)/(2*tI*tI))*(-2+TMath::Erfc((-s+(x[0]-c)*tI/s)/TMath::Sqrt(2)/tI)))); //TMath has an issue that cmath function is used 
    //double f = bl+a/(2*(tD-tI))*tD*tI*((exp((s*s-2*(x[0]-c)*tD)/(2*tD*tD))*(1+erf((-s+(x[0]-c)*tD/s)/sqrt(2)/tD)))+ (exp((s*s-2*(x[0]-c)*tI)/(2*tI*tI))*(-2+erfc((-s+(x[0]-c)*tI/s)/sqrt(2)/tI))));
    return f; 
}

void diff_fit(){
    TCanvas *c1 =new TCanvas("c1","Single signal fit try",1000,800);
    TFile *data = new TFile("thomas.root");
    TH1D *signal = (TH1D *)data->Get("S1_0");
    signal->SetTitle("data 1");
    TF1 *fitF = new TF1("fitF",riseFunc1,0,2900,6);
    fitF->SetParameter(0,0.001);
    fitF->FixParameter(1,22.9765);
    fitF->SetParLimits(2,100000,10000000);
    fitF->SetParLimits(3,1000,2000);
    fitF->SetParLimits(4,0,100);
//  fitF->SetParLimits(5,0,0.0001);
    gStyle->SetOptFit(1011);
    signal->Fit("fitF","R");
}

