#include "CMGTools/MuonCalibration/interface/MuonCalibrator.h"
#include <math.h>


MuonCalibrator::MuonCalibrator(const std::string& magneticFieldMap , const std::string& kalmanCalibration,bool isData) {
  mapFile_ = new TFile(magneticFieldMap.c_str());
  magneticMap_ = (TH2F*)mapFile_->Get("mapCorrection");
  kalmanFile_ = new TFile(kalmanCalibration.c_str());
  kalmanMap_K_ = (TH3F*)kalmanFile_->Get("K_0");
  kalmanMap_A_ = (TH3F*)kalmanFile_->Get("A_0");
  kalmanMap_B_ = (TH3F*)kalmanFile_->Get("B_0");
  kalmanMap_M_ = (TH3F*)kalmanFile_->Get("M_0");
  kalmanCov_ = (TH2F*)kalmanFile_->Get("COVARIANCE_0");
  isData_ = isData;
 }


MuonCalibrator::~MuonCalibrator() {
  mapFile_->Close();
  kalmanFile_->Close();
}


double MuonCalibrator::getCorrectedPt(double pt,double eta,double phi,int charge) {
  double magneticMapFactor = magneticMap_->GetBinContent(magneticMap_->GetBin(
									      magneticMap_->GetXaxis()->FindBin(phi),
									      magneticMap_->GetYaxis()->FindBin(eta)));
  if (!isData_)
    magneticMapFactor=1.0;

  double curvature = magneticMapFactor/pt;
  double sinTheta  = sin(2*atan(exp(-eta))); 


  //  double K = kalmanMap_K_->GetBinContent(kalmanMap_K_->GetBin(                1,
  //									      kalmanMap_K_->GetYaxis()->FindBin(eta),
  //									      kalmanMap_K_->GetZaxis()->FindBin(phi)));
  //  double M = kalmanMap_M_->GetBinContent(kalmanMap_M_->GetBin(                1,
  //									      kalmanMap_M_->GetYaxis()->FindBin(eta),
  //									      kalmanMap_M_->GetZaxis()->FindBin(phi)));


  double K = kalmanMap_K_->GetBinContent(kalmanMap_K_->GetBin(1,1,1));

  double A = kalmanMap_A_->GetBinContent(kalmanMap_A_->GetBin(1,1,1));
  
  double B = kalmanMap_B_->GetBinContent(kalmanMap_B_->GetBin(                1,
									      kalmanMap_B_->GetYaxis()->FindBin(eta),
									      kalmanMap_B_->GetZaxis()->FindBin(phi)));

  double M = kalmanMap_M_->GetBinContent(kalmanMap_M_->GetBin(                1,
									      kalmanMap_M_->GetYaxis()->FindBin(eta),
									      kalmanMap_M_->GetZaxis()->FindBin(phi)));

  curvature = (K*eta*eta+A-1)*curvature +curvature/(1+sinTheta*M*curvature) +charge*B;

  return 1.0/curvature;


}
