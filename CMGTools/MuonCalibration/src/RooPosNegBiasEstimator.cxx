#include "RooFit.h"
#include "Riostream.h"
#include "RooStreamParser.h"
#include "RooNLLVar.h"
#include "RooChi2Var.h"
#include "RooMsgService.h"
#include "TString.h" 
#include "RooRealVar.h"
# include "../interface/RooPosNegBiasEstimator.h"
#include "TMath.h"
#include <math.h>

ClassImp(RooPosNegBiasEstimator)

//_____________________________________________________________________________
  RooPosNegBiasEstimator::RooPosNegBiasEstimator(const char *name, const char *title,RooAbsReal& _bias, int npos,  const Double_t * pos, int nneg, const Double_t * neg):
    RooAbsReal(name,title),
    bias("bias","Observable",this,_bias)
{  
  // Constructor with formula expression and list of input variables


  
  for ( int i=0;i<npos;++i)
    posVector.push_back(pos[i]);
  for ( int i=0;i<nneg;++i)
    negVector.push_back(neg[i]);

  //sort the vectors
  std::sort(posVector.begin(),posVector.end());
  std::sort(negVector.begin(),negVector.end());


  
}





//_____________________________________________________________________________
RooPosNegBiasEstimator::RooPosNegBiasEstimator(const RooPosNegBiasEstimator& other, const char* name) : 
  RooAbsReal(other, name),
  bias("bias",this,other.bias)
{
  posVector = other.posVector;
  negVector = other.negVector;

  // Copy constructor
}



//_____________________________________________________________________________
RooPosNegBiasEstimator::~RooPosNegBiasEstimator() 
{

}







//_____________________________________________________________________________
Double_t RooPosNegBiasEstimator::evaluate() const
{
  return -log(getProbability());


}


Double_t RooPosNegBiasEstimator::getProbability() const
{

  std::vector<double> posNew; 
  std::vector<double> negNew; 


  for (unsigned int i=0;i<posVector.size();++i)
    posNew.push_back(posVector.at(i)+Double_t(bias));
  for (unsigned int j=0;j<negVector.size();++j)
    negNew.push_back(negVector.at(j)-Double_t(bias));

  
  double prob = TMath::KolmogorovTest(posVector.size(),posNew.data(),negVector.size(),negNew.data(),"");
  return prob;


}


//_____________________________________________________________________________
Bool_t RooPosNegBiasEstimator::isValidReal(Double_t /*value*/, Bool_t /*printError*/) const 
{
  // Check if given value is valid
  return kTRUE ;
}



//_____________________________________________________________________________
Bool_t RooPosNegBiasEstimator::redirectServersHook(const RooAbsCollection& newServerList, Bool_t mustReplaceAll, Bool_t nameChange, Bool_t /*isRecursive*/)
{
  // Propagate server change information to embedded RooFormula object
  //  changeDependents(newServerList,mustReplaceAll,nameChange);
  return kTRUE;
}



