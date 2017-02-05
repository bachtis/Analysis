#include "KaMuCa/Derivation/interface/DataSetMaker2D.h"
#include "TLorentzVector.h"
#include "KaMuCa/Calibration/interface/KalmanMuonCalibrator.h"

DataSetMaker2D::DataSetMaker2D(const std::string& outputFileName,const TH3* coords,int bins,double min,double max,int binsy,double miny,double maxy, const char* calib,bool squareSum): 
  coords_(coords),
  squareSum_(squareSum)
{
  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();
  if (strlen(calib)==0) {
    doCalib_=false;
  }
  else {
    calib_ = new KalmanMuonCalibrator(calib);
    doCalib_=true;
  }
  int bin;
  for (int i=1;i<coords->GetNbinsX()+1;++i) {
    for (int j=1;j<coords->GetNbinsY()+1;++j) {
      for (int k=1;k<coords->GetNbinsZ()+1;++k) {
	bin=coords->GetBin(i,j,k);
	char* name;
	if(asprintf(&name,"histo_%d",bin)<0)
	  continue;
	histoMap_[bin] = new TH2D(name,name,bins,min,max,binsy,miny,maxy);
      }
    }
  }


}


void DataSetMaker2D::close(){
  fOut->cd();
  int bin;
  for (int i=1;i<coords_->GetNbinsX()+1;++i) {
    for (int j=1;j<coords_->GetNbinsY()+1;++j) {
      for (int k=1;k<coords_->GetNbinsZ()+1;++k) {
	bin=coords_->GetBin(i,j,k);
	histoMap_[bin]->Write();
      }
    }
  }
  fOut->Close();
}

void DataSetMaker2D::processTree(const std::string&fileName,const std::string& cut,int lepton) {
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
  float eta1;
  float eta2;
  float phi1;
  float phi2;
  float mass;
  float massErr;
  float genMass;
 



  reduced->SetBranchAddress("c1",&c1);
  reduced->SetBranchAddress("c2",&c2);
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

  TLorentzVector v1(1,1,1,1);
  TLorentzVector v2(1,1,1,1);
  int entries=reduced->GetEntries();
  int bin1,bin2,binx,biny,binz;
  double x=0.0;
  double y=0.0;



  for (int i=0;i<entries;++i) {
    reduced->GetEntry(i);
    if (doCalib_) {
      pt1 = calib_->getCorrectedPt(pt1,eta1,phi1,1);
      pt2 = calib_->getCorrectedPt(pt2,eta2,phi2,-1);
    }
    c1=1.0/pt1;
    c2=1.0/pt2;
    v1.SetPtEtaPhiM(pt1,eta1,phi1,0.105);
    v2.SetPtEtaPhiM(pt2,eta2,phi2,0.105);
    mass=(v1+v2).M();
    if (squareSum_)
      binx=coords_->GetXaxis()->FindBin(pt1*pt1+pt2*pt2);
    else
      binx=coords_->GetXaxis()->FindBin(pt1);

    biny=coords_->GetYaxis()->FindBin(eta1);
    binz=coords_->GetZaxis()->FindBin(phi1);

    if (binx==0||binx==coords_->GetNbinsX()+1) 
      continue;
    if (biny==0||biny==coords_->GetNbinsY()+1) 
      continue;
    if (binz==0||binz==coords_->GetNbinsZ()+1) 
      continue;
    

    bin1 = coords_->GetBin(binx,biny,binz);

    if (squareSum_)
      binx=coords_->GetXaxis()->FindBin(pt1*pt1+pt2*pt2);
    else
      binx=coords_->GetXaxis()->FindBin(pt2);
    biny=coords_->GetYaxis()->FindBin(eta2);
    binz=coords_->GetZaxis()->FindBin(phi2);

    if (binx==0||binx==coords_->GetNbinsX()+1) 
      continue;
    if (biny==0||biny==coords_->GetNbinsY()+1) 
      continue;
    if (binz==0||binz==coords_->GetNbinsZ()+1) 
      continue;

    bin2 = coords_->GetBin(binx,biny,binz);


    x=mass;
    y=massErr;
    if (lepton==1)
      histoMap_[bin1]->Fill(x,y);

    if (lepton==2)
      histoMap_[bin2]->Fill(x,y);

    if (lepton==3 && bin1==bin2)
      histoMap_[bin2]->Fill(x,y);


    if (i % 1000000==0)
	printf("Processed %d \%d entries\n",i,entries);

  }

  fOut->cd();
  fIn->Close();

}




