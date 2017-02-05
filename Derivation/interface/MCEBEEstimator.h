#ifndef KaMuCa_Derivation_MCEBEEstimator_h
#define KaMuCa_Derivation_MCEBEEstimator_h


#include "TFile.h"
#include "TTree.h"
#include "TH3D.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
class KalmanMuonCalibrator;
class MCEBEEstimator  {
 public:
  MCEBEEstimator(const std::string& outputFileName,const TH3D* map,const char* calib = "MC_Moriond17_13TeV");
  void processTree(const std::string& fileName,const std::string& cut);
  void close();
 private:
  KalmanMuonCalibrator *calib_;
  TH3D* resolution_;
  TH3D* resolutionRef_;
  TFile *fOut;
};



#endif
