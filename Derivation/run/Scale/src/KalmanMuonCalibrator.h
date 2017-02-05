#ifndef KalmanMuonCalibrator_h
#define KalmanMuonCalibrator_h


#include "SimpleKalmanCalculator.h"
#include "TH1D.h"
#include "TH2D.h"
#include <math.h>
#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
class KalmanMuonCalibrator {
 public:
  KalmanMuonCalibrator();
  ~KalmanMuonCalibrator();

  void processFile(const char*,const char*);
  void processFileMC(const char*,const char*);
  void save(const char*);
  void load(const char*);
  enum  Measurement {A,K,L,A11,A12,A21,A22,A31,A32,e,e2,B0,B11,B12,B21,B22,B31,B32};


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
  void processLine(double c1,double eta1,double phi1,double c2,double eta2,double phi2,double scale,double resolution,int Z);
  SimpleKalmanCalculator * calculator_;
  TMatrixDSym correlationMatrix();

  
};


#endif
