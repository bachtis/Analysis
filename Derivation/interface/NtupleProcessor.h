#ifndef KaMuCa_Derivation_NtupleProcessor_h
#define KaMuCa_Derivation_NtupleProcessor_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

class NtupleProcessor  {
 public:
  NtupleProcessor(const std::string& outputFileName,bool isData,float target0=3.095,float target1=0.0,float width=0.0);
  void processTree(const std::string& fileName,const std::string& cut);
  void write();
 private:
  RooWorkspace *w;
  RooDataSet *data;
  TFile  *fOut;
  bool isData;
  float width_;
  float target0_;
  float target1_;
  KalmanMuonCalibrator *calib_;
};



#endif
