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
#include "TH1F.h"

ClassImp(RooPosNegBiasEstimator)

//_____________________________________________________________________________
  RooPosNegBiasEstimator::RooPosNegBiasEstimator(const char *name, const char *title,RooAbsReal& _bias, const RooDataSet* pos, const RooDataSet* neg,const char* var1, const char* var2):
    RooAbsReal(name,title),
    bias("bias","Observable",this,_bias)
{  
  // Constructor with formula expression and list of input variables

  for (int i=0;i<pos->numEntries();++i) {
    posVector.push_back(((RooRealVar*)pos->get(i)->find(var1))->getVal());
    posWeights.push_back(pos->weight());
  }
  for (int i=0;i<neg->numEntries();++i) {
    negVector.push_back(((RooRealVar*)neg->get(i)->find(var2))->getVal());
    negWeights.push_back(neg->weight());
  }
  
}





//_____________________________________________________________________________
RooPosNegBiasEstimator::RooPosNegBiasEstimator(const RooPosNegBiasEstimator& other, const char* name) : 
  RooAbsReal(other, name),
  bias("bias",this,other.bias)
{
  posVector = other.posVector;
  negVector = other.negVector;
  posWeights = other.posVector;
  negWeights = other.negVector;
  
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

  TH1F * hPos = new TH1F("pos","pos",25,0.005,0.05);
  hPos->Sumw2();
  TH1F * hNeg = new TH1F("neg","neg",25,0.005,0.05);
  hNeg->Sumw2();

  for (unsigned int i=0;i<posVector.size();++i) {
    hPos->Fill(posVector[i]+Double_t(bias),posWeights[i]);
  }
  for (unsigned int j=0;j<negVector.size();++j) {
    hNeg->Fill(negVector[j]-Double_t(bias),negWeights[j]);
  }
  
  //  double prob = TMath::KolmogorovTest(posVector.size(),posNew.data(),negVector.size(),negNew.data(),"");
  double prob = hPos->Chi2Test(hNeg,"WW");
  delete hPos;
  delete hNeg;
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



