
#include "TProfile3D.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "TH2D.h"
#include "TFile.h" 
#include "TRandom2.h"

class UnfoldingMap {
 public:
  UnfoldingMap(int binspt,double minpt,double maxpt,int binssigma,double minsigma,double maxsigma);
  void loadSpectra(const std::string& spectraFile, const std::string& spectrum); 
  void process(int N1,int N2);
  const TProfile3D * getMap() const {
    return profile;
  }
  void write(); 
 

  ~UnfoldingMap();

 private:
  TProfile3D * profile;
  TFile *fOut_;
  TFile *fIn_;
  TH2D *spectrum_;
  RooDataHist *datahist;
  RooHistPdf *histpdf;
  RooRealVar *x;
  RooRealVar *y;
  TRandom *random_;
  double errMin_;
  double errMax_;
};
