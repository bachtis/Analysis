#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CMGTools/MuonCalibration/interface/UnfoldingMap.h"
#include <math.h>
#include "TString.h"


UnfoldingMap::UnfoldingMap(int binscurv,double mincurv,double maxcurv,int binssigma,double minsigma,double maxsigma)  {

  fOut_ = new TFile("mapResult.root","RECREATE");
  fOut_->cd();
  profile = new TProfile3D("unfolding","unfolding",binscurv,mincurv,maxcurv,binscurv,mincurv,maxcurv,binssigma,minsigma,maxsigma);
  x = new RooRealVar("x","x",1/10.);
  y = new RooRealVar("y","y",1/10.);
  random_ = new TRandom2(1010821982);
  errMin_ = minsigma;
  errMax_ = maxsigma;
 }



void UnfoldingMap::loadSpectra(const std::string& spectraFile, const std::string& spectrum) {
  edm::FileInPath path(spectraFile);

  fIn_ = new TFile(path.fullPath().c_str());
  spectrum_ = (TH2D*)fIn_->Get(spectrum.c_str());
  datahist = new RooDataHist("hist","hist",RooArgList(*x,*y),spectrum_);
  histpdf  = new RooHistPdf("histpdf","hist",RooArgSet(*x,*y),*datahist,2);
} 


UnfoldingMap::~UnfoldingMap() {
  fIn_->Close();
    if(datahist)
      delete datahist;
    if(histpdf)
      delete histpdf;
    if(x)
      delete x;
    if(y)
      delete y;
    if(random_)
      delete random_;
}



void UnfoldingMap::write() {
  fOut_->cd();
  profile->Write();
  fOut_->Close();
}

void UnfoldingMap::process(int N1,int N2) {
  int samples = N1/100000;
  double xv=0.0;
  double yv=0.0;
  double sxv=0.0;
  double syv=0.0;
  double errx=0.0;
  double erry=0.0;
  double errMx=0.0;
  double errMy=0.0;
  for( int i=0;i<samples;++i) {
    printf("Starting next 100000 events\n");
    RooDataSet *data = histpdf->generate(RooArgSet(*x,*y),100000);
    for ( int j=0;j<data->numEntries();++j) {
      const RooArgSet* line = data->get(j);
      xv = ((RooRealVar*)line->find("x"))->getVal();
      yv = ((RooRealVar*)line->find("y"))->getVal();
      for ( int k=0;k<N2;++k) {
	errx = (errMin_+random_->Rndm()*(errMax_-errMin_));
	errMx =errx*xv; 
	sxv = xv+random_->Gaus(0.0,errMx);
	erry = (errMin_+random_->Rndm()*(errMax_-errMin_))*yv;
	errMy =erry*yv; 
	syv = yv+random_->Gaus(0.0,errMy);
	profile->Fill(sxv,syv,errx,xv-sxv);
	profile->Fill(syv,sxv,erry,yv-syv);

      }

    }
    if (data)
      delete data;
    
  }


}


