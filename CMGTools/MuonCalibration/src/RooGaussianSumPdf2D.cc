#include "../interface/RooGaussianSumPdf2D.h"
#include "TMath.h"
#include "RooMath.h"
#include "DataFormats/Math/interface/LorentzVector.h"

ClassImp(RooGaussianSumPdf2D)

RooGaussianSumPdf2D::	RooGaussianSumPdf2D(const char *name, const char *title,
			      RooAbsReal& _u1,
			      RooAbsReal& _u2,
			      RooAbsReal& _scale1,
			      RooAbsReal& _scale2,
			      RooAbsReal& _error1,
			      RooAbsReal& _error2,
		              const RooDataSet& _data,
			      const char* varName1,
			      const char* varName2):

RooAbsPdf(name,title),
u1("u1","u1",this,_u1),
u2("u2","u2",this,_u2),
scale1("scale1","scale1",this,_scale1),
scale2("scale2","scale2",this,_scale2),
error1("error1","error1",this,_error1),
error2("error2","error2",this,_error2)
{
  double weight = 1.0;
  for(  int i=0;i<_data.numEntries();++i) {
    const RooArgSet * line = _data.get(i);
    if (_data.isWeighted())
      weight = _data.weight();
    weights.push_back(weight);
    data1.push_back(((RooAbsReal*)line->find(varName1))->getVal());
    data2.push_back(((RooAbsReal*)line->find(varName2))->getVal());
  }


}


RooGaussianSumPdf2D::RooGaussianSumPdf2D(const RooGaussianSumPdf2D& other, const char* name) :
RooAbsPdf(other,name),
u1("u1",this,other.u1),
u2("u2",this,other.u2),
scale1("scale1",this,other.scale1),
scale2("scale2",this,other.scale2),
error1("error1",this,other.error1),
error2("error2",this,other.error2),
data1(other.data1),
data2(other.data2),
weights(other.weights)
{
 

}


Double_t RooGaussianSumPdf2D::evaluate() const
{
  Double_t arg1= 0.0;
  Double_t arg2= 0.0;
  Double_t sum=0.0;
  Double_t sumw=0.0;
  static const Double_t sqrt2pi(sqrt(2*TMath::Pi()));

  for (unsigned int i=0;i<data1.size();++i) {
    arg1 = (u1-scale1*data1.at(i))/error1;
    arg2 = (u2-scale2*data2.at(i))/error2;
    sum=sum+(weights.at(i)/(sqrt2pi*error1*error2))*exp(-0.5*arg1*arg1-0.5*arg2*arg2);
    sumw=sumw+weights.at(i);
  }
  return sum/sumw;
}


