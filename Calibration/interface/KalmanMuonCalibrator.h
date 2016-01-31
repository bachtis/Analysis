#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1I.h"
#include "TProfile2D.h"
#include "TRandom3.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
class KalmanMuonCalibrator {
 public:
  KalmanMuonCalibrator();
  KalmanMuonCalibrator(const std::string&);
  double getCorrectedPt(double pt,double eta,double phi,int charge);
  double smear(double pt,double eta);
  double smearUsingEbE(double pt,double eta,double error);
  double smearForSync(double pt,double eta);
  double getCorrectedPtMag(double,double,double);
  double getCorrectedError(double pt,double eta,double error);
  double getCorrectedErrorAfterSmearing(double pt,double eta,double error);

  int getN();
  void vary(int,int);
  void varyClosure(int);
  void reset();
  void randomize();


  ~KalmanMuonCalibrator();



 private:
  double closure(double,double);
  TRandom * random_;
  void resetHisto(TH1*,const TH1* );
  int varyClosure_;

  bool isData_;
  TFile *file_;
  TH2F *magnetic; 


  TH3F *scale_A1; 
  TH3F *scale_A2; 
  TH3F *scale_e;
  TH3F *scale_B0;
  TH3F *scale_B1;
  TH3F *scale_B2;
  TH3F *scale_C1;
  TH3F *scale_C2;

  TH3F *shifted_A1; 
  TH3F *shifted_A2; 
  TH3F *shifted_e;
  TH3F *shifted_B0;
  TH3F *shifted_B1;
  TH3F *shifted_B2;
  TH3F *shifted_C1;
  TH3F *shifted_C2;


  TH1D* aSRC_;
  TH1D* bSRC_;
  TH1D* cSRC_;
  TH1D* dSRC_;
  TH1D* aTARGET_;
  TH1D* bTARGET_;
  TH1D* cTARGET_;
  TH1D* dTARGET_;

  TH1D* aEBE_;
  TH1D* bEBE_;
  TH1D* cEBE_;
  TH1D* dEBE_;




  TH3F *closure_;
  TMatrixD *eigenvectors_;
  TVectorD *eigenvalues_;
  TH1I *covHistoMap_;
  TH1I *covBinMap_;
};
