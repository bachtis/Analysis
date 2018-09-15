#include "KaMuCa/Derivation/interface/NtupleProcessor.h"
#include "TLorentzVector.h"
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"
#include "TRandom.h"
NtupleProcessor::  NtupleProcessor(const std::string& outputFileName,bool isdata,float target0,float target1,float width,const char* calib,int fullCalib) {
  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();
  calib_ = new KalmanMuonCalibrator(calib);
  isData=isdata;
  width_ = width;    
  target0_ = target0;
  target1_ = target1;
  fullCalib_ = fullCalib;
}


void NtupleProcessor::close(){
  fOut->Close();
}

void NtupleProcessor::processTree(const std::string&fileName,const std::string& cut,double maxRap) {
  TFile *fIn = new TFile(fileName.c_str()); 
  TTree *t = (TTree*)fIn->Get("tree");
  //reduce it!
  fOut->cd();
  TTree* reduced = t->CopyTree(cut.c_str());
  reduced->SetName("reduced");

  float pt1;
  float pt2;
  float c1;
  float c2;
  float gc1;
  float gc2;

  float gcSmeared1;
  float gcSmeared2;

  float cErr1;
  float cErr2;
  
  float eta1;
  float eta2;
  float phi1;
  float phi2;
  float mass;
  float massErr;
  float genMass;
  float scale;
  float resolution;
  float rapidity;
  float genMassSmeared;


  reduced->SetBranchAddress("c1",&c1);
  reduced->SetBranchAddress("c2",&c2);
  reduced->SetBranchAddress("gc1",&gc1);
  reduced->SetBranchAddress("gc2",&gc2);

  reduced->SetBranchAddress("cErr1",&cErr1);
  reduced->SetBranchAddress("cErr2",&cErr2);
  reduced->SetBranchAddress("pt1",&pt1);
  reduced->SetBranchAddress("pt2",&pt2);
  reduced->SetBranchAddress("eta1",&eta1);
  reduced->SetBranchAddress("eta2",&eta2);
  reduced->SetBranchAddress("phi1",&phi1);
  reduced->SetBranchAddress("phi2",&phi2);
  reduced->SetBranchAddress("mass",&mass);
  reduced->SetBranchAddress("genMass",&genMass);
  reduced->SetBranchAddress("massErr",&massErr);

  data = reduced->CloneTree(0);
  data->SetName("tree");

  data->Branch("scale",&scale,"scale/F");

  data->Branch("resolution",&resolution,"resolution/F");
  data->Branch("rapidity",&rapidity,"rapidity/F");

  if (!isData) {
    data->Branch("genMassSmeared",&genMassSmeared,"genMassSmeared/F");
    data->Branch("gcSmeared1",&gcSmeared1,"gcSmeared1/F");
    data->Branch("gcSmeared2",&gcSmeared2,"gcSmeared2/F");
  }

  TLorentzVector v1(1,1,1,1);
  TLorentzVector v2(1,1,1,1);
  int entries=reduced->GetEntries();
  double trueMass=0.0;
  
  TRandom *random = new TRandom(10101982);

  for (int i=0;i<entries;++i) {
    reduced->GetEntry(i);
    v1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);
    rapidity =(v1+v2).Rapidity();
    if (fabs(rapidity)>maxRap)
      continue;

    if (fullCalib_==2) {
      pt1 = calib_->getCorrectedPt(pt1,eta1,phi1,1);
      pt2 = calib_->getCorrectedPt(pt2,eta2,phi2,-1);
      cErr1 = calib_->getCorrectedError(pt1,eta1,cErr1);
      cErr2 = calib_->getCorrectedError(pt2,eta2,cErr2);

    }
    else if (fullCalib_==1) {
      pt1 = calib_->getCorrectedPt(pt1,eta1,phi1,1);
      pt2 = calib_->getCorrectedPt(pt2,eta2,phi2,-1);
    }
    else {
      pt1 = calib_->getCorrectedPtMag(pt1,eta1,phi1);
      pt2 = calib_->getCorrectedPtMag(pt2,eta2,phi2);
    }
    c1=1.0/pt1;
    c2=1.0/pt2;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);
    mass=(v1+v2).M();
    
    if (fullCalib_==2) {
      massErr = 0.5*mass*sqrt(cErr1*cErr1+cErr2*cErr2);
    }


    trueMass = target0_+target1_*rapidity*rapidity*rapidity*rapidity;



    if (!isData) {
      genMassSmeared = genMass+random->Gaus(0.0,massErr);
      gcSmeared1 = gc1+random->Gaus(0.0,cErr1*gc1);
      gcSmeared2 = gc2+random->Gaus(0.0,cErr2*gc2);
    }


    scale=mass/trueMass;
    resolution = sqrt((massErr*massErr)/(trueMass*trueMass)+width_*width_/(trueMass*trueMass));

    //    scale=genMassSmeared/trueMass;


    data->Fill();

    if (i % 1000000==0)
	printf("Processed %d \%d entries\n",i,entries);

  }

  fOut->cd();
  data->Write();
  fIn->Close();

}




