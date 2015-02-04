
#include "TProfile3D.h"
#include "TProfile2D.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "RooDerivative.h"
#include "TH2D.h"
#include "TFile.h" 
#include "TRandom2.h"

class UnfoldingMapDeriv {
 public:
  UnfoldingMapDeriv(int binsDeriv,double minDeriv,double maxDeriv,int binssigma,double minsigma,double maxsigma);
  void loadSpectra(const std::string& spectraFile, const std::string& spectrum); 
  void process(int N1,int N2);
  const TProfile2D * getMap() const {
    return profile;
  }
  void write(); 
 
  ~UnfoldingMapDeriv();

 private:
  TProfile2D * profile;
  TFile *fOut_;
  TFile *fIn_;
  TH2D *spectrum_;
  RooDataHist *datahist;
  RooHistPdf *histpdf;
  RooDerivative * dx_;
  RooDerivative * dy_;
  RooRealVar *x;
  RooRealVar *y;
  TRandom *random_;
  double errMin_;
  double errMax_;
};
