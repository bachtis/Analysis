#include "../interface/RooGaussianSumPdfWithSigma2.h"
#include "TMath.h"
#include "RooMath.h"
#include "DataFormats/Math/interface/LorentzVector.h"

ClassImp(RooGaussianSumPdfWithSigma2)

RooGaussianSumPdfWithSigma2::	RooGaussianSumPdfWithSigma2(const char *name, const char *title,
			      RooAbsReal& _mass,
			      RooAbsReal& _scale,
			      RooAbsReal& _error1,
			      RooAbsReal& _error2,
			      const RooDataSet& _data,
			      const char* varName,
 			      const char* errorVarName,
                              const char* errorVarName2): 
RooAbsPdf(name,title),
mass("mass","mass",this,_mass),
scale("scale","scale",this,_scale),
error1("error1","error1",this,_error1),
error2("error2","error2",this,_error2)
{
  double weight = 1.0;
  for(  int i=0;i<_data.numEntries();++i) {
    const RooArgSet * line = _data.get(i);
    if (_data.isWeighted())
      weight = _data.weight();
    weights.push_back(weight);
    masses.push_back(((RooAbsReal*)line->find(varName))->getVal());
    errors1.push_back(((RooAbsReal*)line->find(errorVarName))->getVal());
    errors2.push_back(((RooAbsReal*)line->find(errorVarName2))->getVal());
  }


}


RooGaussianSumPdfWithSigma2::RooGaussianSumPdfWithSigma2(const RooGaussianSumPdfWithSigma2& other, const char* name) :
RooAbsPdf(other,name),
mass("mass",this,other.mass),
scale("scale",this,other.scale),
error1("error1",this,other.error1),
error2("error2",this,other.error2),
weights(other.weights),
masses(other.masses),
errors1(other.errors1),
errors2(other.errors2)
{
 

}


Double_t RooGaussianSumPdfWithSigma2::evaluate() const
{
  Double_t arg= 0.0;
  Double_t sum=0.0;
  Double_t sumw=0.0;
  
  for (unsigned int i=0;i<masses.size();++i) {
    arg = (mass-scale*masses[i])/sqrt(error1*errors1[i]*error1*errors1[i]+error2*errors2[i]*error2*errors2[i]);
    sum=sum+weights[i]*exp(-0.5*arg*arg);
    sumw=sumw+weights[i];
  }
  return sum/sumw;
}

Int_t RooGaussianSumPdfWithSigma2::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,mass)) return 1 ;
  return 0;
}


Double_t RooGaussianSumPdfWithSigma2::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  static const Double_t root2 = sqrt(2.) ;
  static const Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
  Double_t ret = 0;
  Double_t sumw = 0;
  Double_t xscale = 0.0;

  for (unsigned int i=0;i<weights.size();++i) {
    xscale = root2*sqrt(error1*errors1[i]*error1*errors1[i]+error2*errors2[i]*error2*errors2[i]);

    ret = ret+weights[i]*rootPiBy2*sqrt(error1*errors1[i]*error1*errors1[i]+error2*errors2[i]*error2*errors2[i])*(RooMath::erf((mass.max(rangeName)-scale*masses[i])/xscale)-RooMath::erf((mass.min(rangeName)-scale*masses[i])/xscale));
    sumw=sumw+weights[i];
  }
  return ret/sumw ;

}

