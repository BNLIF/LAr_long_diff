//fit rising edge with simple error function

//fit parameters:
//                  par0 = rise edge amplitude;
//                  par1 = drift time;
//                  par2 = rise time;
//                  par3 = base line shift;
double riseFunc2(double *x, double *par){
	double a = par[0];// edge amplitude
        double c = par[1];// drift time
	double s = par[2];// rise time
	double bl = par[3];//base line shift
	double f = bl + a*(TMath::Erf((x[0]-c)/TMath::Sqrt(2)/s)+1)/2;
	return f; 
}        
void diff_erf(){
    TCanvas *c1 =new TCanvas("c1","Single signal fit try",1000,800);
    TFile *data = new TFile("thomas.root");
    TH1D *signal = (TH1D *)data->Get("S1_0");
    signal->SetTitle("data 1");
    TF1 *fitF = new TF1("fitF",riseFunc2,1400,2200,4);
    fitF->SetParameter(0,0.001);
    //fitF->FixParameter(1,22.9765);
    fitF->SetParLimits(1,1500,2000);
    fitF->SetParLimits(2,0,20000);
    fitF->SetParLimits(3,0,1);
    gStyle->SetOptFit(1111);
    signal->Fit("fitF","R");
}

