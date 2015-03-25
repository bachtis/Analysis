#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooAbsCategory.h"
#include "RooDataSet.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "TH1.h"
#include "RooGaussian.h"
#include "TRandom.h"




class DynamicBinnedSmearingPdf : public RooAbsPdf {
public:
	DynamicBinnedSmearingPdf() {} ;
	DynamicBinnedSmearingPdf(const char *name, const char *title,
				 RooAbsReal& _mass,
				 RooAbsReal& _scale,
				 RooAbsReal& _error,
				 const TH1F* histo);
	DynamicBinnedSmearingPdf(const DynamicBinnedSmearingPdf& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new DynamicBinnedSmearingPdf(*this,newname); }
	inline virtual ~DynamicBinnedSmearingPdf() { }

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
	Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

protected:
	
	RooRealProxy mass;
	RooRealProxy scale;
	RooRealProxy error;
	Double_t evaluate() const ;
	TH1F *data_;
	Double_t resolution(double x) const  {
	  return exp(-0.5*(x-scale)*(x-scale)/(error*error))/(sqrt(2*M_PI)*error); 
	}



private:
	ClassDef(DynamicBinnedSmearingPdf,1) // Your description goes here...                                                                                                   
	  
};
