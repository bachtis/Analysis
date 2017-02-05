#include "KaMuCa/Derivation/interface/MCEBEEstimator.h"
#include "TLorentzVector.h"
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"
#include "TRandom.h"
MCEBEEstimator::  MCEBEEstimator(const std::string& outputFileName,const TH3D* hmap,const char* calib) {
  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();
  calib_ = new KalmanMuonCalibrator(calib);
  resolution_ = new TH3D(*hmap);
  resolution_->SetName("resolution");
  resolutionRef_ = new TH3D(*hmap);
  resolutionRef_->SetName("resolutionRef");
  calib_ = new KalmanMuonCalibrator(calib);
}


void MCEBEEstimator::close(){
  fOut->cd();
  resolution_->Write();
  resolutionRef_->Write();
  fOut->Close();
}

void MCEBEEstimator::processTree(const std::string&fileName,const std::string& cut) {
  TFile *fIn = new TFile(fileName.c_str()); 
  TTree *t = (TTree*)fIn->Get("tree");
  //reduce it!
  fOut->cd();
  TTree* reduced = t->CopyTree(cut.c_str());
  reduced->SetName("reduced");

  float pt1;
  float pt2;
  float c1;
  float c2;
  float cErr1;
  float cErr2;

  float gc1;
  float gc2;
  float eta1;
  float eta2;
  float phi1;
  float phi2;

  float gpt1,gpt2;
  
  reduced->SetBranchAddress("c1",&c1);
  reduced->SetBranchAddress("c2",&c2);
  reduced->SetBranchAddress("gc1",&gc1);
  reduced->SetBranchAddress("gc2",&gc2);
  reduced->SetBranchAddress("pt1",&pt1);
  reduced->SetBranchAddress("pt2",&pt2);
  reduced->SetBranchAddress("eta1",&eta1);
  reduced->SetBranchAddress("eta2",&eta2);
  reduced->SetBranchAddress("phi1",&phi1);
  reduced->SetBranchAddress("phi2",&phi2);
  reduced->SetBranchAddress("cErr1",&cErr1);
  reduced->SetBranchAddress("cErr2",&cErr2);   
  TLorentzVector v1(1,1,1,1);
  TLorentzVector v2(1,1,1,1);
  int entries=reduced->GetEntries();
  TRandom *random = new TRandom(10101982+resolution_->GetEntries());
  for (int i=0;i<entries;++i) {
    reduced->GetEntry(i);
    pt1 = calib_->getCorrectedPt(pt1,eta1,phi1,1);
    pt2 = calib_->getCorrectedPt(pt2,eta2,phi2,-1);
    gpt1=1.0/gc1;
    gpt2=1.0/gc2;
    c1=1.0/pt1;
    c2=1.0/pt2;

    resolution_->Fill(gpt1,eta1,c1/gc1);
    resolution_->Fill(gpt2,eta2,c2/gc2);

    //now the smeared
    c1=random->Gaus(gc1,cErr1*gc1);
    c2=random->Gaus(gc2,cErr2*gc2);
 
    resolutionRef_->Fill(gpt1,eta1,c1/gc1);
    resolutionRef_->Fill(gpt2,eta2,c2/gc2);



    if (i % 1000000==0)
	printf("Processed %d \%d entries\n",i,entries);

  }
  fIn->Close();

}




