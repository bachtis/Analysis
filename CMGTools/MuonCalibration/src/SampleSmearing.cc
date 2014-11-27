#include "CMGTools/MuonCalibration/interface/SampleSmearing.h"
#include <math.h>
#include "TString.h"
#include "RooRealVar.h"
SampleSmearing::SampleSmearing(const std::string& inputFile) {
  inputFile_ = new TFile(inputFile.c_str());
  inputData_ = (RooDataSet*)inputFile_->Get("data");
  outputFile_ = new TFile("large_"+TString(inputFile.c_str()),"RECREATE");
  random_ = new TRandom3(1010821982);
 }


SampleSmearing::~SampleSmearing() {
  inputFile_->Close();
}


void SampleSmearing::smear(int N) {

  double c1=0.0;
  double c2=0.0;
  double mErr1=0.0;
  double mErr2=0.0;
  double mass=0.0;

  outputFile_->cd();
  RooDataSet * data2 = new RooDataSet("data","data",*inputData_->get());

  for (int i=0; i < inputData_->numEntries();++i) {
      const RooArgSet*  line = inputData_->get(i);

    c1 = ((RooRealVar*)line->find("curvRaw1"))->getVal();
    c2 = ((RooRealVar*)line->find("curvRaw2"))->getVal();
    mErr1 = ((RooRealVar*)line->find("massErrRaw1"))->getVal();
    mErr2 = ((RooRealVar*)line->find("massErrRaw2"))->getVal();
    mass = ((RooRealVar*)line->find("massRaw"))->getVal();
    
    for (int j=0;j<N;++j) {
      c1=c1+random_->Gaus(0.0,2*mErr1*c1/mass);
      c2=c2+random_->Gaus(0.0,2*mErr2*c2/mass);
      ((RooRealVar*)line->find("curvRaw1"))->setVal(c1);
      ((RooRealVar*)line->find("curvRaw2"))->setVal(c2);
      data2->add(*line);
    }
  }

  outputFile_->cd();
  data2->Write();
  outputFile_->Close();
}
