#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "RooDataSet.h"
 

class SampleSmearing {
 public:
  SampleSmearing(const std::string& file);
  void smear(int N);
  ~SampleSmearing();


 private:
  TFile *inputFile_;
  TFile *outputFile_;
  RooDataSet * inputData_; 
  TRandom *random_;

};
