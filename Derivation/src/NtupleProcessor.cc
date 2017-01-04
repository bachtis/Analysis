#include "KaMuCa/Derivation/interface/NtupleProcessor.h"
#include "TLorentzVector.h"

NtupleProcessor::NtupleProcessor(const std::string& outputFileName,bool isdata,float target0,float target1,float width) {

  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();

  isData=isdata;
  target0_=target0;
  target1_=target1;
  width_=width;
  

  w= new RooWorkspace("w","w");
  w->factory("c1[0,1]");
  w->factory("c2[0,1]");
  w->factory("eta1[-2.6,2.6]");
  w->factory("eta2[-2.6,2.6]");
  w->factory("phi1[-4,4]");
  w->factory("phi2[-4,4]");
  w->factory("scale[0,100000]");
  w->factory("resolution[0,100]");

  w->defineSet("vars","c1,c2,eta1,eta2,phi1,phi2,scale,resolution");
  data=new RooDataSet("data","data",*w->set("vars"));

  if(isData)
    printf("Will Calibrate the magnetic field using previous calibrations\n");

  calib_ = new KalmanMuonCalibrator("DATA_80X_13TeV");
}



void NtupleProcessor::processTree(const std::string&fileName,const std::string& cut) {
  TFile *fIn = new TFile(fileName.c_str()); 
  TTree *t = (TTree*)fIn->Get("tree");
  //reduce it!
  fOut->cd();
  TTree* reduced = t->CopyTree(cut.c_str());
  
  float pt1;
  float pt2;
  float c1;
  float c2;
  float eta1;
  float eta2;
  float phi1;
  float phi2;
  float mass;
  float massErr;


  reduced->SetBranchAddress("c1",&c1);
  reduced->SetBranchAddress("c2",&c2);
  reduced->SetBranchAddress("pt1",&pt1);
  reduced->SetBranchAddress("pt2",&pt2);
  reduced->SetBranchAddress("eta1",&eta1);
  reduced->SetBranchAddress("eta2",&eta2);
  reduced->SetBranchAddress("phi1",&phi1);
  reduced->SetBranchAddress("phi2",&phi2);
  reduced->SetBranchAddress("mass",&mass);
  reduced->SetBranchAddress("massErr",&massErr);

  TLorentzVector v1(1,1,1,1);
  TLorentzVector v2(1,1,1,1);
  float scale=0.0;
  float resolution=0.0;
  float trueMass=0.0;
  int entries=reduced->GetEntries();
  for (int i=0;i<entries;++i) {
    reduced->GetEntry(i);
    w->var("eta1")->setVal(eta1);
    w->var("eta2")->setVal(eta2);
    w->var("phi1")->setVal(phi1);
    w->var("phi2")->setVal(phi2);
    if (isData) {
      pt1 = calib_->getCorrectedPtMag(pt1,eta1,phi1);
      pt2 = calib_->getCorrectedPtMag(pt2,eta2,phi2);
    }
      c1=1.0/pt1;
      c2=1.0/pt2;
      w->var("c1")->setVal(c1);
      w->var("c2")->setVal(c2);
      v1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
      v2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);
      mass=(v1+v2).M();
      trueMass=target0_+target1_*(v1+v2).Rapidity()*(v1+v2).Rapidity();
      scale=mass/trueMass;
      resolution =sqrt(massErr*massErr+width_*width_)/mass;
      w->var("scale")->setVal(scale);
      w->var("resolution")->setVal(resolution);
      RooArgSet set(*w->set("vars"));
      data->add(set);

      if (i % 1000000==0)
	printf("Processed %d \%d entries\n",i,entries);

  }
  delete reduced;
  fIn->Close();

}



void NtupleProcessor::write() {
  fOut->cd();
  data->Write();
  fOut->Close();
}
