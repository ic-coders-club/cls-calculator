#include <iostream>
#include <vector>

#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"

//Jad Marrouche, November 2013
//to compile, do:
//g++ CLs_calculator.cpp -I`root-config --incdir` `root-config --libs` -o CLs_calculator.exe

int main(void) {

  //int nbins=3;
  int nbins=1;

  std::vector<double> bg_orig(nbins, 0.0);
  std::vector<double> bguncert_orig(nbins, 0.0);
  std::vector<double> sig_orig(nbins, 0.0);
  std::vector<double> siguncert_orig(nbins, 0.0);
  std::vector<double> data(nbins, 0.0);

  //THINGS TO CHANGE START
  //bg_orig[0] = 1.1; bg_orig[1] = 1.2; bg_orig[2] = 2.6;
  //bguncert_orig[0] = 1.0; bguncert_orig[1] = 1.0; bguncert_orig[2] = 0.54;
  //sig_orig[0] = 2.44; sig_orig[1] = 2.96; sig_orig[2] = 4.81;
  //siguncert_orig[0] = 0.2; siguncert_orig[1] = 0.2; siguncert_orig[2] = 0.20;
  //data[0] = 1; data[1] = 0; data[2] = 3;
  bg_orig[0]=50.0;
  bguncert_orig[0]=0.10;
  sig_orig[0]=10.0;
  siguncert_orig[0]=0.2;
  data[0]=60;
  //THINGS TO CHANGE END


  TH1D * histbg = new TH1D("histbg", "", 1000, -49.95, 49.05);
  TH1D * histsig = new TH1D("histsig", "", 1000, -49.95, 49.05);
  
  TH1D * histsigbg = new TH1D("histsigbg", "", 1000, -49.95, 49.05);

  TH1D * CLs = new TH1D("CLs","", 1000, -49.95, 49.05);

  TRandom3 * randomnumber = new TRandom3();
  
  int num_iterations = 10000;

  std::vector<double> bglambda(nbins, 0.0);
  std::vector<double> siglambda(nbins, 0.0);
  std::vector<int> bg(nbins, 0);
  std::vector<int> sigbg(nbins, 0);
  
  for(int i=0; i<num_iterations; i++) {
  
    double qsigbg=1.0, qbg=1.0;
    
    for(int j=0; j<nbins; j++) {
      bglambda[j] = randomnumber->TRandom::Gaus(bg_orig[j], bguncert_orig[j] * bg_orig[j]);
      siglambda[j] = randomnumber->TRandom::Gaus(sig_orig[j], siguncert_orig[j] * sig_orig[j]);
      bg[j] = randomnumber->TRandom::Poisson(bglambda[j]);
      sigbg[j] = randomnumber->TRandom::Poisson(bglambda[j] + siglambda[j]);
      qbg *= TMath::Poisson(bg[j], bg_orig[j]+sig_orig[j]) / TMath::Poisson(bg[j], bg_orig[j]);
      qsigbg *= TMath::Poisson(sigbg[j], bg_orig[j]+sig_orig[j]) / TMath::Poisson(sigbg[j], bg_orig[j]);
    }
    
    histbg->Fill(-2.0 * log(qbg));
    histsigbg->Fill(-2.0 * log(qsigbg));
		 
  }

  double qbgd=1.0;

  for(int j=0; j<nbins; j++) {
    qbgd *= TMath::Poisson(data[j], bg_orig[j]+sig_orig[j]) / TMath::Poisson(data[j], bg_orig[j]);
  }
  
  qbgd = -2.0 * log(qbgd);

  std::cout << "-2*ln(Q) = " << qbgd << std::endl;


  double runningsumbg=histbg->Integral();
  double runningsumsigbg=histsigbg->Integral();
  
  for(int i=0; i<CLs->GetNbinsX(); i++) {
    if(runningsumbg > 0.0) {
      CLs->SetBinContent(i+1, runningsumsigbg / runningsumbg);
    }
    runningsumbg -= histbg->GetBinContent(i+1);
    runningsumsigbg -= histsigbg->GetBinContent(i+1);
  }

  std::cout << "CLs = " << CLs->GetBinContent(CLs->FindBin(qbgd)) << std::endl;

  TFile * ofile = new TFile("output-file.root","RECREATE");
  ofile->cd();
  histbg->Write();
  histsigbg->Write();
  CLs->Write();
  ofile->Close();

  return 0;
}
