#include "../interface/RooGaussianSumPdfWithSigma.h"
#include "TMath.h"
#include "RooMath.h"
#include "DataFormats/Math/interface/LorentzVector.h"

ClassImp(RooGaussianSumPdfWithSigma)

RooGaussianSumPdfWithSigma::	RooGaussianSumPdfWithSigma(const char *name, const char *title,
			      RooAbsReal& _mass,
			      RooAbsReal& _scale,
			      RooAbsReal& _error1,
			      RooAbsReal& _error2,
			      const RooDataSet& _data,
			      const char* varName,
							   const char* errorVarName): 
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
    errors.push_back(((RooAbsReal*)line->find(errorVarName))->getVal());
  }


}


RooGaussianSumPdfWithSigma::RooGaussianSumPdfWithSigma(const RooGaussianSumPdfWithSigma& other, const char* name) :
RooAbsPdf(other,name),
mass("mass",this,other.mass),
scale("scale",this,other.scale),
error1("error1",this,other.error1),
error2("error2",this,other.error2),
weights(other.weights),
masses(other.masses),
errors(other.errors)
{
 

}


Double_t RooGaussianSumPdfWithSigma::evaluate() const
{
  Double_t arg= 0.0;
  Double_t error= 0.0;
  Double_t sum=0.0;
  Double_t sumw=0.0;
  static const Double_t rootPiTimes2 = sqrt(2*atan2(0.0,-1.0));
  
  
  for (unsigned int i=0;i<masses.size();++i) {
    error = errors[i]*errors[i]+error1*(masses[i]*masses[i])+error2;
    //protection!
    if (error<0.0)
      continue;
    else
      error=sqrt(error);
    arg = (mass-scale*masses[i])/(error);
    sum=sum+weights[i]*exp(-0.5*arg*arg)/(rootPiTimes2*error);
    sumw=sumw+weights[i];
  }
  return sum/sumw;
}

Int_t RooGaussianSumPdfWithSigma::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,mass)) return 1 ;
  return 0;
}


Double_t RooGaussianSumPdfWithSigma::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  static const Double_t root2 = sqrt(2.) ;
  //  static const Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
  Double_t ret = 0;
  Double_t sumw = 0;
  Double_t xscale = 0.0;
  Double_t error = 0.0;
  
  for (unsigned int i=0;i<weights.size();++i) {
    error = (errors[i]*errors[i]+error1*(masses[i]*masses[i])+error2);
    if (error<0.0)
      continue;
    else
      error=sqrt(error);

    xscale = root2*(error);

    ret = ret+weights[i]*0.5*(RooMath::erf((mass.max(rangeName)-scale*masses[i])/xscale)-RooMath::erf((mass.min(rangeName)-scale*masses[i])/xscale));
    sumw=sumw+weights[i];
  }
  return ret/sumw ;

}

