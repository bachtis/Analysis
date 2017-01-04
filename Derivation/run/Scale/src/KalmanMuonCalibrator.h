#ifndef KalmanMuonCalibrator_h
#define KalmanMuonCalibrator_h


#include "SimpleKalmanCalculator.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "RooArgSet.h"
#include "TAxis.h"
#include "TFile.h"
#include "RooDataSet.h"
class KalmanMuonCalibrator {
 public:
  KalmanMuonCalibrator();
  ~KalmanMuonCalibrator();

  void processFile(const char*,const char*);
  enum  Measurement {A,K,L,A11,A12,A21,A22,A31,A32,e,B0,B11,B12,B21,B22,B31,B32};

 private:

  double *data_;
  double *derivative_;
  unsigned int nEtaBins_;
  TAxis* etaAxis_;
  unsigned int nMaterialEtaBins_;
  TAxis* etaMaterialAxis_;

  unsigned int totalBins_;

  unsigned int getBin(Measurement measurement,float eta);
  double getData(Measurement measurement,float eta);
  void setData(Measurement measurement,float eta,double data);
  void addDerivative(Measurement measurement,float eta,double data);
  void resetDerivative();
  void processLine(const RooArgSet*);
  SimpleKalmanCalculator * calculator_;


  
};


#endif
