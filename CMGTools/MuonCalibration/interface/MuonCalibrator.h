#include "TFile.h"
#include "TH3F.h"
#include "TH2F.h"


class MuonCalibrator {
 public:
  MuonCalibrator(const std::string& magneticFieldMap , const std::string& kalmanCalibration,bool isData = false);
  double getCorrectedPt(double,double,double,int);

  ~MuonCalibrator();



 private:
  bool isData_;
  TFile *mapFile_;
  TFile *kalmanFile_;
  TH2F *magneticMap_; 
  TH3F *kalmanMap_K_; 
  TH3F *kalmanMap_A_; 
  TH3F *kalmanMap_B_; 
  TH3F *kalmanMap_M_; 
  TH2F  *kalmanCov_;
};
