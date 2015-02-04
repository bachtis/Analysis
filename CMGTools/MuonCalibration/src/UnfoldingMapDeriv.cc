#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CMGTools/MuonCalibration/interface/UnfoldingMapDeriv.h"
#include <math.h>
#include "TString.h"


UnfoldingMapDeriv::UnfoldingMapDeriv(int binsDeriv,double minDeriv,double maxDeriv,int binssigma,double minsigma,double maxsigma)  {

  fOut_ = new TFile("mapResult.root","RECREATE");
  fOut_->cd();
  profile = new TProfile2D("unfolding","unfolding",binsDeriv,minDeriv,maxDeriv,binssigma,minsigma,maxsigma);
  x = new RooRealVar("x","x",1/10.);
  y = new RooRealVar("y","y",1/10.);
  random_ = new TRandom2(1010821982);
  errMin_ = minsigma;
  errMax_ = maxsigma;
 }



void UnfoldingMapDeriv::loadSpectra(const std::string& spectraFile, const std::string& spectrum) {
  edm::FileInPath path(spectraFile);

  fIn_ = new TFile(path.fullPath().c_str());
  spectrum_ = (TH2D*)fIn_->Get(spectrum.c_str());
  datahist = new RooDataHist("hist","hist",RooArgList(*x,*y),spectrum_);
  histpdf  = new RooHistPdf("histpdf","hist",RooArgSet(*x,*y),*datahist,2);
  dx_ = histpdf->derivative(*x,1,0.001);
  dy_ = histpdf->derivative(*y,1,0.001);
} 


UnfoldingMapDeriv::~UnfoldingMapDeriv() {
  fIn_->Close();
    if(datahist)
      delete datahist;
    if(histpdf)
      delete histpdf;
    if(dx_)
      delete dx_;
    if(dy_)
      delete dy_;
    if(x)
      delete x;
    if(y)
      delete y;
    if(random_)
      delete random_;
}



void UnfoldingMapDeriv::write() {
  fOut_->cd();
  profile->Write();
  fOut_->Close();
}

void UnfoldingMapDeriv::process(int N1,int N2) {

  int samples = N1/100000;
  double xv=0.0;
  double yv=0.0;
  double sxv=0.0;
  double syv=0.0;
  double dx=0.0;
  double dy=0.0;
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

        //calculate the derivatives:
        x->setVal(sxv);
        y->setVal(syv);
        dx = dx_->getVal()/histpdf->getVal();
        dy = dy_->getVal()/histpdf->getVal();

	profile->Fill(dx,errx,xv-sxv);
	profile->Fill(dy,erry,yv-syv);
      }
      
    }
    if (data)
      delete data;
    
 }


}


