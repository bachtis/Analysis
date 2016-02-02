#ifndef KaMuCa_Derivation_NtupleProcessorLuca_h
#define KaMuCa_Derivation_NtupleProcessorLuca_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"

class NtupleProcessorLuca  {
 public:
  NtupleProcessorLuca(const std::string& outputFileName,bool isData,bool isGen);
  void processTree(const std::string& fileName,const std::string& cut);
  void write();
 private:
  RooWorkspace *w;
  RooDataSet *data;
  TFile  *fOut;
  bool isData;
  bool isGen;

  

};



#endif
