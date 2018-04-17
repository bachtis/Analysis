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
#include "RooDataHist.h"
#include "TH1F.h"
#include "TRandom.h"




class RooGaussianSumPdfRelative : public RooAbsPdf {
public:
	RooGaussianSumPdfRelative() {} ;
	RooGaussianSumPdfRelative(const char *name, const char *title,
		      RooAbsReal& _mass,
		      RooAbsReal& _scale,
		      RooAbsReal& _error1,
		      RooAbsReal& _error2,
 		      const RooAbsData& data,
		      const char* varName);



	RooGaussianSumPdfRelative(const RooGaussianSumPdfRelative& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooGaussianSumPdfRelative(*this,newname); }
	inline virtual ~RooGaussianSumPdfRelative() { }

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
	Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;


protected:
	
	RooRealProxy mass ;
	RooRealProxy scale ;
	RooRealProxy error1 ;
	RooRealProxy error2 ;
	
	Double_t evaluate() const ;
	std::vector<std::pair<double,double> > data;


private:


	ClassDef(RooGaussianSumPdfRelative,1) // Your description goes here...                                                                                                   
	  
};
