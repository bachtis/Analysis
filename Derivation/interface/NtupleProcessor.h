#ifndef KaMuCa_Derivation_NtupleProcessor_h
#define KaMuCa_Derivation_NtupleProcessor_h


#include "TFile.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"

class NtupleProcessor  {
 public:
  NtupleProcessor(const std::string& outputFileName,bool isData,bool isGen,const std::string& treePrefix = "Onia");

  void processTree(const std::string& fileName,const std::string& cut);
  void write();
 private:
  RooWorkspace *w;
  RooDataSet *data;
  TFile  *fOut;
  bool isData;
  bool isGen;
  std::string prefix;
  

};



#endif
