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




class RooGaussianSumPdf2D : public RooAbsPdf {
public:
	RooGaussianSumPdf2D() {} ;
	RooGaussianSumPdf2D(const char *name, const char *title,
		      RooAbsReal& _u1,
		      RooAbsReal& _u2,
		      RooAbsReal& _scale1,
		      RooAbsReal& _scale2,
		      RooAbsReal& _error1,
		      RooAbsReal& _error2,
 		      const RooDataSet& data,
  	              const char* varName1,
		      const char* varName2
			  );



	RooGaussianSumPdf2D(const RooGaussianSumPdf2D& other, const char* name=0) ;
	virtual TObject* clone(const char* newname) const { return new RooGaussianSumPdf2D(*this,newname); }
	inline virtual ~RooGaussianSumPdf2D() { }

	//	Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
	//	Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;

	Bool_t selfNormalized() const {
	  return true;
	}
protected:
	
	RooRealProxy u1 ;
	RooRealProxy u2 ;
	RooRealProxy scale1 ;
	RooRealProxy scale2 ;
	RooRealProxy error1 ;
	RooRealProxy error2 ;
	
	Double_t evaluate() const ;
	std::vector<double> data1;
	std::vector<double> data2;
        std::vector<double> weights;


private:


	ClassDef(RooGaussianSumPdf2D,1) // Your description goes here...                                                                                                   
	  
};
