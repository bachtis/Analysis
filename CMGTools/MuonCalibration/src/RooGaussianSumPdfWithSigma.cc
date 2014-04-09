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
			      const char* sigmaName):

RooAbsPdf(name,title),
mass("mass","mass",this,_mass),
scale("scale","scale",this,_scale),
error1("error1","error1",this,_error1),
error2("error2","error2",this,_error2)
{
  for(  int i=0;i<_data.numEntries();++i) {
    const RooArgSet * line = _data.get(i);
    if (_data.isWeighted())
      weight.push_back(_data.weight());
    else
      weight.push_back(1.);
    data.push_back(((RooAbsReal*)line->find(varName))->getVal());
    sigma.push_back(((RooAbsReal*)line->find(sigmaName))->getVal());
  }


}


RooGaussianSumPdfWithSigma::RooGaussianSumPdfWithSigma(const RooGaussianSumPdfWithSigma& other, const char* name) :
RooAbsPdf(other,name),
mass("mass",this,other.mass),
scale("scale",this,other.scale),
error1("error1",this,other.error1),
error2("error2",this,other.error2),
data(other.data),
sigma(other.sigma),
weight(other.weight)

{
 

}


Double_t RooGaussianSumPdfWithSigma::evaluate() const
{
  Double_t arg= 0.0;
  Double_t sum=0.0;
  Double_t sumw=0.0;
  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));

  for (unsigned int i=0;i<data.size();++i) {
    arg = (mass-scale*data.at(i))/(error1*sigma.at(i));
    sum=sum+(weight.at(i)/(sqrt2pi*error1*sigma.at(i)))*exp(-0.5*arg*arg);
    sumw=sumw+weight.at(i);
  }
  return sum/sumw;
}


