#include <TROOT.h>
#include "TApplication.h"
#include "Rtypes.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "math.h"
#include "stdio.h"
#include "TH1F.h"
#include "TString.h"
#include <vector>
#include "TChain.h"
#include "TFile.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCut.h"
#include "TLeaf.h"
#include "TH2F.h"
#include "TMath.h"

using namespace::std;
int main (Int_t argc, char *argv[])
{
  ifstream infile("thomas.data");
  TString filename,dir;
  Double_t dis;
  Double_t E;
  Int_t fileno = 0;

  TString temp;
  
  TFile *file = new TFile("thomas.root","RECREATE");
  TTree *T = new TTree("T","T");
  T->SetDirectory(file);
  T->Branch("fileno",&fileno,"fileno/I");
  T->Branch("E",&E,"E/D");
  T->Branch("dis",&dis,"dis/D");

  while (!infile.eof()){
    infile >> dir >> filename >> dis;
    if (dis == -1) break;
    TString input = "./processed/thomas/" + filename;
    Double_t time[35000],data[16][35000];
    Double_t HV[20];
    Int_t ncount = 0 ;
    Int_t nline;
    cout << input << endl;
    ifstream infile1(input);
    infile1 >> nline >> dis;
    cout << nline << " " << dis << endl;
    if (nline!=1){
      infile1 >> temp;
      for (Int_t i=0;i!=nline;i++){
   	infile1 >> HV[i];
        cout << HV[i]<< endl;
      }
    }else{
      infile1 >> temp >> temp;
    }
    
    while (!infile1.eof()){
      infile1 >> time[ncount];
      time[ncount]*=1000.;
      // if (ncount!=0){
      // 	time[ncount] = time[ncount-1] +=1.25; 
      // }
      for (Int_t i=0;i!=nline;i++){
  	infile1 >> data[i][ncount];
      }
      if (data[0][ncount]!=0 || time[ncount]!=0){
	ncount++;
      }
    }
    ncount --;

    for (Int_t i=0;i!=nline;i++){
      TH1D *h1 = new TH1D(Form("S1_%d",fileno),Form("S1_%d",fileno),ncount/2.,time[0],time[ncount-1]);
      TH1D *h4 = new TH1D(Form("S4_%d",fileno),Form("S4_%d",fileno),ncount/2.,time[0],time[ncount-1]);
      for (Int_t j=0;j!=ncount;j++){
  	h1->Fill(time[j],data[i][j]);
  	h4->Fill(time[j],1);
      }
      h1->Divide(h4);
      h1->SetDirectory(file);
      h4=0;
      E = HV[i];
      T->Fill();
      fileno ++;
    }
    infile1.close();
  }
  infile.close();
  file->Write();
  file->Close();

}
