#ifndef KaMuCa_Derivation_DataSetMaker2D_h
#define KaMuCa_Derivation_DataSetMaker2D_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "TH3.h"
#include "TH1D.h"

class KalmanMuonCalibrator;
class DataSetMaker2D  {
 public:
  DataSetMaker2D(const std::string& outputFileName,const TH3* coords, int binsx,double minx,double maxx,int binsy,double miny,double maxy, const char* calib = "DATA_Moriond17_13TeV",bool squareSum = false);
  void processTree(const std::string& fileName,const std::string& cut,int lepton = 1);
  void close();
 private:
  TTree *data;
  TFile  *fOut;
  KalmanMuonCalibrator *calib_;
  const TH3* coords_;
  std::map<int,TH2D*> histoMap_;
  bool doCalib_;
  bool squareSum_;
};



#endif
