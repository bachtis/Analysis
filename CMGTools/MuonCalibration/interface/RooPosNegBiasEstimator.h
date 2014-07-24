#ifndef ROO_POSNEGBIASESTIMATOR
#define ROO_kPOSNEGBIASESTIMATOR


#include "RooAbsReal.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooLinkedList.h"
#include "RooAbsProxy.h"
#include <TMath.h>
#include "RooDataSet.h"
class RooArgSet ;

class RooPosNegBiasEstimator : public RooAbsReal {
public:
  // Constructors, assignment etc
  inline RooPosNegBiasEstimator()   { }
  RooPosNegBiasEstimator(const char *name, const char *title,RooAbsReal& _bias, const RooDataSet* pos, const RooDataSet* neg,const char* var1,const char* var2);
    RooPosNegBiasEstimator(const RooPosNegBiasEstimator& other, const char* name=0);
    virtual TObject* clone(const char* newname) const { return new RooPosNegBiasEstimator(*this,newname); }
    virtual ~RooPosNegBiasEstimator();
    Double_t getProbability() const;


protected:

  // Function evaluation
  virtual Double_t evaluate() const ;
  // Post-processing of server redirection
  virtual Bool_t redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t isRecursive) ;

  virtual Bool_t isValidReal(Double_t value, Bool_t printError) const ;
  RooRealProxy bias;
  std::vector<double> posVector;
  std::vector<double> posWeights;
  std::vector<double> negVector;
  std::vector<double> negWeights;
  ClassDef(RooPosNegBiasEstimator,1) // Real-valued function of other RooAbsArgs calculated by a TFormula expression
};

#endif
