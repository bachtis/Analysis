#include "../interface/RooGaussianSumPdfRelative.h"
#include "TMath.h"
#include "RooMath.h"
#include "DataFormats/Math/interface/LorentzVector.h"

ClassImp(RooGaussianSumPdfRelative)

RooGaussianSumPdfRelative::	RooGaussianSumPdfRelative(const char *name, const char *title,
			      RooAbsReal& _mass,
			      RooAbsReal& _scale,
			      RooAbsReal& _error1,
			      RooAbsReal& _error2,
			      const RooAbsData& _data,
 			      const char* varName):

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
    std::pair<double,double> pair = std::make_pair(weight,((RooAbsReal*)line->find(varName))->getVal());
    data.push_back(pair);
  }


}


RooGaussianSumPdfRelative::RooGaussianSumPdfRelative(const RooGaussianSumPdfRelative& other, const char* name) :
RooAbsPdf(other,name),
mass("mass",this,other.mass),
scale("scale",this,other.scale),
error1("error1",this,other.error1),
error2("error2",this,other.error2),
data(other.data)

{
 

}


Double_t RooGaussianSumPdfRelative::evaluate() const
{
  Double_t arg= 0.0;
  Double_t sum=0.0;
  Double_t sumw=0.0;
  
  for (unsigned int i=0;i<data.size();++i) {
    arg = (mass-scale*data[i].second)/(error1*data[i].second);
    sum=sum+data[i].first*exp(-0.5*arg*arg);
    sumw=sumw+data[i].first;
  }
  return sum/sumw;
}

Int_t RooGaussianSumPdfRelative::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  if (matchArgs(allVars,analVars,mass)) return 1 ;
  return 0;
}


Double_t RooGaussianSumPdfRelative::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  static const Double_t root2 = sqrt(2.) ;
  static const Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
  Double_t xscale = root2*error1;
  Double_t ret = 0;
  Double_t sumw = 0;
  
  for (unsigned int i=0;i<data.size();++i) {
    ret = ret+data[i].first*rootPiBy2*(error1*data[i].second)*(RooMath::erf((mass.max(rangeName)-scale*data[i].second)/(xscale*data[i].second))-RooMath::erf((mass.min(rangeName)-scale*data[i].second)/(xscale*data[i].second)));
    sumw=sumw+data[i].first;
  }
  return ret/sumw ;

}

