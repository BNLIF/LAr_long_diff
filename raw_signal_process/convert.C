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
int
main (Int_t argc, char *argv[])
{
//void convert(){
  ifstream infile("output.txt");
  TString dir1, dir2, shaping;
  Double_t temp;
  Int_t fileno;
  Double_t HV, dis;
  TFile *file;
  TTree *T;
  while (!infile.eof()){
    infile >> dir1 >> dir2 >> temp;
    TString rootfile;
    rootfile = dir1 + "_" + dir2 + ".root";
    file = new TFile(rootfile,"RECREATE");
    T = new TTree("T","T");
    T->SetDirectory(file);
    T->Branch("fileno",&fileno,"fileno/I");
    T->Branch("HV",&HV,"HV/D");
    T->Branch("dis",&dis,"dis/D");

    Double_t time[20000],data[3][20000];
    Int_t ncount = 0 ;
    if (!dir1.IsNull()){
      TString filename;
      filename = "./processed/" + dir1 + "/" + dir2 + "/list.dat";
      ifstream infile1(filename);
      cout << filename << endl;
      while(!infile1.eof()){
      	infile1 >> fileno >> shaping >> HV >> dis;
	if (!shaping.IsNull()){
	  T->Fill();
	  filename = "./processed/" + dir1 + "/" + dir2 + Form("/%d.dat",fileno);
	  ifstream infile2(filename);
	  infile2 >> shaping >> HV >> dis;
	  ncount = 0;
	  while(!infile2.eof()){
	    infile2 >> time[ncount] >> data[0][ncount] >> data[1][ncount] >> data[2][ncount];
	    time[ncount]*=1e9;
	    ncount++;
	  }
	  ncount --;
	  
	  TH1D *h1 = new TH1D(Form("S1_%d",fileno),Form("S1_%d",fileno),ncount,time[0],time[ncount-1]);
	  TH1D *h2 = new TH1D(Form("S2_%d",fileno),Form("S2_%d",fileno),ncount,time[0],time[ncount-1]);
	  TH1D *h3 = new TH1D(Form("S3_%d",fileno),Form("S3_%d",fileno),ncount,time[0],time[ncount-1]);
	  TH1D *h4 = new TH1D(Form("S4_%d",fileno),Form("S4_%d",fileno),ncount,time[0],time[ncount-1]);
	  for (Int_t i=0;i!=ncount;i++){
	    h1->Fill(time[i],data[0][i]);
	    h2->Fill(time[i],data[1][i]);
	    h3->Fill(time[i],data[2][i]);
	    h4->Fill(time[i],1);
	  }
	  h1->Divide(h4);
	  h2->Divide(h4);
	  h3->Divide(h4);
	  h1->SetDirectory(file);
	  h2->SetDirectory(file);
	  h3->SetDirectory(file);

	  
	  
	}
      }
      infile1.close();
    }
    file->Write();
    file->Close();
    //break;
  }
}
