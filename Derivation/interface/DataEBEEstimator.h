#ifndef KaMuCa_Derivation_DataEBEEstimator_h
#define KaMuCa_Derivation_DataEBEEstimator_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TH3.h"
#include "TH1D.h"

class KalmanMuonCalibrator;
class DataEBEEstimator  {
 public:
  DataEBEEstimator(const std::string& outputFileName,const TH3* coords,int genMass,int bins,double min,double max, const char* calib = "DATA_Moriond17_13TeV");
  void processTree(const std::string& fileName,const std::string& cut);
  void close();
 private:
  TTree *data;
  TFile  *fOut;
  KalmanMuonCalibrator *calib_;
  const TH3* coords_;
  std::map<int,TH1D*> histoMap_;
  bool doCalib_;
  int doGenMass_;

};



#endif
