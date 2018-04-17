#include "KaMuCa/Derivation/interface/DATAEBEEstimator.h"
#include "TLorentzVector.h"
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

DATAEBEEstimator::DATAEBEEstimator(const std::string& outputFileName,const TH3* coords,int binsMass,double minMass,double maxMass,int binsPt,double minPt,double maxPt, const char* calib): 
  coords_(coords)
{
  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();
  if (strlen(calib)==0) {
    doCalib_=false;
  }
  else {
    calib_ = new KalmanMuonCalibrator(calib);
    doCalib_=true;
  }

  for (int i=1;i<coords->GetNbinsX()+1;++i) {
    for (int j=1;j<coords->GetNbinsY()+1;++j) {
      for (int k=1;k<coords->GetNbinsZ()+1;++k) {
	int bin=coords->GetBin(i,j,k);
	char* name;
	if(asprintf(&name,"histo_%d",bin)<0)
	  continue;
	histoMap_[bin] = new TH3D(name,name,binsMass,minMass,maxMass,binsPt,minPt,maxPt,binsPt,minPt,maxPt);
      }
    }
  }


}



void DATAEBEEstimator::close(){
  fOut->cd();
  for (int i=1;i<coords_->GetNbinsX()+1;++i) {
    for (int j=1;j<coords_->GetNbinsY()+1;++j) {
      for (int k=1;k<coords_->GetNbinsZ()+1;++k) {
	int bin=coords_->GetBin(i,j,k);
	histoMap_[bin]->Write();
      }
    }
  }


  fOut->Close();
}

void DATAEBEEstimator::processTree(const std::string&fileName) {
  TFile *fIn =new TFile(fileName.c_str()); 
  TTree *reduced = (TTree*)fIn->Get("tree");
  
  float pt1;
  float pt2;
  float c1;
  float c2;
  float eta1;
  float eta2;
  float phi1;
  float phi2;
  float mass;
  float massErr;


  reduced->SetBranchAddress("c1",&c1);
  reduced->SetBranchAddress("c2",&c2);
  reduced->SetBranchAddress("pt1",&pt1);
  reduced->SetBranchAddress("pt2",&pt2);
  reduced->SetBranchAddress("eta1",&eta1);
  reduced->SetBranchAddress("eta2",&eta2);
  reduced->SetBranchAddress("phi1",&phi1);
  reduced->SetBranchAddress("phi2",&phi2);
  reduced->SetBranchAddress("mass",&mass);
  reduced->SetBranchAddress("massErr",&massErr);

  TLorentzVector v1(1,1,1,1);
  TLorentzVector v2(1,1,1,1);
  int entries=reduced->GetEntries();
  int bin1,bin2;

  int binx,biny,binz;
  for (int i=0;i<entries;++i) {
    reduced->GetEntry(i);
    if (doCalib_) {
      pt1 = calib_->getCorrectedPt(pt1,eta1,phi1,1);
      pt2 = calib_->getCorrectedPt(pt2,eta2,phi2,-1);
    }

    c1=1.0/pt1;
    c2=1.0/pt2;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);
    mass=(v1+v2).M();

    binx=coords_->GetXaxis()->FindBin(pt1);
    biny=coords_->GetYaxis()->FindBin(eta1);
    binz=coords_->GetZaxis()->FindBin(phi1);

    if (binx==0 || binx>coords_->GetNbinsX()||biny==0 || biny>coords_->GetNbinsY()  ||binz==0 || binz>coords_->GetNbinsZ())
      continue;
  
    bin1=coords_->GetBin(binx,biny,binz);

    binx=coords_->GetXaxis()->FindBin(pt2);
    biny=coords_->GetYaxis()->FindBin(eta2);
    binz=coords_->GetZaxis()->FindBin(phi2);
    if (binx==0 || binx>coords_->GetNbinsX()||biny==0 || biny>coords_->GetNbinsY()  ||binz==0 || binz>coords_->GetNbinsZ())
      continue;

    bin2=coords_->GetBin(binx,biny,binz);


    if (bin1 !=bin2)
      continue;

    histoMap_[bin1]->Fill(mass,pt1,pt2);

    if (i % 1000000==0)
	printf("Processed %d \%d entries\n",i,entries);

  }

  fOut->cd();
  fIn->Close();

}




