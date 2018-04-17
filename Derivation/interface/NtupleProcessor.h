#ifndef KaMuCa_Derivation_NtupleProcessor_h
#define KaMuCa_Derivation_NtupleProcessor_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
class KalmanMuonCalibrator;
class NtupleProcessor  {
 public:
  NtupleProcessor(const std::string& outputFileName,bool isData,float target0=3.095,float target1=0.0,float width=0.0,const char* calib = "DATA_Moriond17_13TeV",int fullCalib = 0);
  void processTree(const std::string& fileName,const std::string& cut,double maxRap = 5.0);
  void close();
 private:
  RooWorkspace *w;
  TTree *data;
  TFile  *fOut;
  bool isData;
  float width_;
  float target0_;
  float target1_;
  KalmanMuonCalibrator *calib_;
  int fullCalib_;
};



#endif
