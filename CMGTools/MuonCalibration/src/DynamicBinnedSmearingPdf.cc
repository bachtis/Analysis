#include "../interface/DynamicBinnedSmearingPdf.h"
#include "TMath.h"
#include "RooMath.h"
#include "DataFormats/Math/interface/LorentzVector.h"

ClassImp(DynamicBinnedSmearingPdf)

DynamicBinnedSmearingPdf::DynamicBinnedSmearingPdf(const char *name, const char *title,
			      RooAbsReal& _mass,
			      RooAbsReal& _scale,
			      RooAbsReal& _error,
                              const TH1F* data):

RooAbsPdf(name,title),
mass("mass","mass",this,_mass),
scale("scale","scale",this,_scale),
error("error","error",this,_error)
{
  data_ = (TH1F*)data->Clone();
  data_->SetName("internalH");
  data_->Scale(1./data_->Integral());
}


DynamicBinnedSmearingPdf::DynamicBinnedSmearingPdf(const DynamicBinnedSmearingPdf& other, const char* name) :
RooAbsPdf(other,name),
mass("mass",this,other.mass),
scale("scale",this,other.scale),
error("error",this,other.error),
data_(other.data_)
{
 

}


Double_t DynamicBinnedSmearingPdf::evaluate() const
{
  int bin = data_->GetXaxis()->FindBin(mass);
  if (bin<1 || bin >data_->GetNbinsX())
    return 0.0;
  double sum=0.0;
  for (int i=1;i<data_->GetNbinsX()+1;++i) {
    sum=sum+resolution(data_->GetXaxis()->GetBinCenter(bin)-data_->GetXaxis()->GetBinCenter(i))*data_->GetBinContent(i);
  }


  return sum;

}

Int_t DynamicBinnedSmearingPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* /*rangeName*/) const 
{
  //  if (matchArgs(allVars,analVars,mass)) return 1 ;
  //  return 0;
  return 1;
}



Double_t DynamicBinnedSmearingPdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
    double sum=0.0;
    //    int binmin=data_->GetXaxis()->FindBin(mass.min(rangeName));
    //    int binmax=data_->GetXaxis()->FindBin(mass.max(rangeName));

    for (int j=1; j<data_->GetNbinsX()+1; ++j) 
      for (int i=1;i<data_->GetNbinsX()+1;++i) 
	sum=sum+resolution(data_->GetXaxis()->GetBinCenter(j)-data_->GetXaxis()->GetBinCenter(i))*data_->GetBinContent(i)*data_->GetXaxis()->GetBinWidth(i);
    return sum;
}

