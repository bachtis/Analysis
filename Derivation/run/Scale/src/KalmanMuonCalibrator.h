#ifndef KalmanMuonCalibrator_h
#define KalmanMuonCalibrator_h
#include "TMatrixDSymfwd.h"
class SimpleKalmanCalculator;
class TAxis;


class KalmanMuonCalibrator {
 public:
  KalmanMuonCalibrator();
  ~KalmanMuonCalibrator();

  void processFile(const char*,const char*);
  void processFileMC(const char*,const char*);
  void save(const char*);
  void load(const char*);


  enum  Measurement {A,K,e,B,L};


 private:
  SimpleKalmanCalculator * calculator_;

  double *data_;
  double *derivative_;

  unsigned int nEtaMaterialBins_;
  TAxis* etaMaterialAxis_;

  unsigned int nEtaBins_;
  TAxis* etaAxis_;

  unsigned int nPhiBins_;
  TAxis* phiAxis_;

  unsigned int totalBins_;



  unsigned int getBin(Measurement measurement,float eta,float phi);
  double getData(Measurement measurement,float eta,float phi);
  void setData(Measurement measurement,float eta,float phi,double data);
  void addDerivative(Measurement measurement,float eta,float phi,double data);
  void resetDerivative();
  void processLine(double c1,double eta1,double phi1,double c2,double eta2,double phi2,double scale,double resolution,int Z);
  TMatrixDSym correlationMatrix();


  
};


#endif
