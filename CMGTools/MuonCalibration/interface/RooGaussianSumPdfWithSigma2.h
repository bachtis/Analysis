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
#include "RooDataSet.h"
#include "TH1F.h"
#include "TRandom.h"




class RooGaussianSumPdfWithSigma2 : public RooAbsPdf {
public:
	RooGaussianSumPdfWithSigma2() {} ;
	RooGaussianSumPdfWithSigma2(const char *name, const char *title,
		      RooAbsReal& _mass,
		      RooAbsReal& _scale,
		      RooAbsReal& _error1,
		      RooAbsReal& _error2,
 		      const RooDataSet& data,
		      const char* varName,
 	              const char* errorVarName,
		      const char* errorVarName2);



	RooGaussianSumPdfWithSigma2(const RooGaussianSumPdfWithSigma2& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooGaussianSumPdfWithSigma2(*this,newname); }
	inline virtual ~RooGaussianSumPdfWithSigma2() { }

	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
	Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;


protected:
	
	RooRealProxy mass ;
	RooRealProxy scale ;
	RooRealProxy error1 ;
	RooRealProxy error2 ;
	
	Double_t evaluate() const ;
	std::vector<double> weights;
	std::vector<double> masses;
	std::vector<double> errors1;
	std::vector<double> errors2;



private:


	ClassDef(RooGaussianSumPdfWithSigma2,1) // Your description goes here...                                                                                                   
	  
};
