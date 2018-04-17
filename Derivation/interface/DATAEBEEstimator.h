#ifndef KaMuCa_Derivation_DataEBEEstimator_h
#define KaMuCa_Derivation_DataEBEEstimator_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TH3.h"
#include "TH3D.h"

class KalmanMuonCalibrator;
class DATAEBEEstimator  {
 public:
  DATAEBEEstimator(const std::string& outputFileName,const TH3* eta, int binsMass,double minMass,double maxMass,int binsPT,double minPT,double maxPT, const char* calib = "DATA_Moriond17_13TeV");
  void processTree(const std::string& fileName);
  void close();
 private:
  TTree *data;
  TFile  *fOut;
  KalmanMuonCalibrator *calib_;
  const TH3* coords_;
  std::map<int,TH3D*> histoMap_;
  bool doCalib_;
};



#endif
