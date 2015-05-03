void merge(){
  TFile *file = new TFile("result_thomas_new.root");
  
  TTree *T = (TTree*)file->Get("T");
  Int_t fileno;
  Double_t E, dis, diff_time, drift_time;
  Double_t cutoff, m_charge;
  Double_t m_diff_time, m_diff_time_err;
  Double_t m_drift_time;
  
  T->SetBranchAddress("fileno",&fileno);
  T->SetBranchAddress("E",&E);
  T->SetBranchAddress("dis",&dis);
  T->SetBranchAddress("diff_time",&diff_time);
  T->SetBranchAddress("drift_time",&drift_time);
  T->SetBranchAddress("cutoff",&cutoff);
  T->SetBranchAddress("m_charge",&m_charge);
  T->SetBranchAddress("m_diff_time",&m_diff_time);
  T->SetBranchAddress("m_diff_time_err",&m_diff_time_err);
  T->SetBranchAddress("m_drift_time",&m_drift_time);

  ifstream infile("./results/result.txt");
  Double_t range_min[200],range_max[200];
  Double_t fit_diff[200], fit_diff_err[200];
  Double_t chi2ndf[200],vel[200],vel_err[200];
  Int_t temp;
  for (Int_t i=0;i!=185;i++){
    infile >> temp;
    infile >> range_min[temp] >> range_max[temp] 
	   >> fit_diff[temp] >> fit_diff_err[temp]
	   >> chi2ndf[temp] >> vel[temp] >> vel_err[temp];
    vel[temp] -= 65;
  }
  
  
  TFile *file1 = new TFile("merge.root","RECREATE");
  TTree *t1 = new TTree("T","T");
  t1->SetDirectory(file1);
  t1->Branch("fileno",&fileno,"fileno/I");
  t1->Branch("E",&E,"E/D");
  t1->Branch("dis",&dis,"dis/D");
  t1->Branch("diff_time",&diff_time,"diff_time/D");
  t1->Branch("drift_time",&drift_time,"drift_time/D");
  t1->Branch("cutoff",&cutoff,"cutoff/D");
  t1->Branch("m_charge",&m_charge,"m_charge/D");
  t1->Branch("m_diff_time",&m_diff_time,"m_diff_time/D");
  t1->Branch("m_diff_time_err",&m_diff_time_err,"m_diff_time_err/D");
  t1->Branch("m_drift_time",&m_drift_time,"m_drift_ime/D");

  Double_t m_drift_time_err;
  t1->Branch("m_drift_time_err",&m_drift_time_err,"m_drift_ime_err/D");
  
  
  Double_t x2ndf;
  Double_t min,max;
  t1->Branch("x2ndf",&x2ndf,"x2ndf/D");
  t1->Branch("min",&min,"min/D");
  t1->Branch("max",&max,"max/D");
  
  for (Int_t i=0;i!=T->GetEntries();i++){
    T->GetEntry(i);
    m_diff_time = fit_diff[fileno];
    m_diff_time_err = fit_diff_err[fileno];
    m_drift_time = vel[fileno];
    m_drift_time_err = vel_err[fileno];
    min = range_min[fileno];
    max = range_max[fileno];
    x2ndf = chi2ndf[fileno];
    t1->Fill();
  }
  file1->Write();
  file1->Close();
  

}
