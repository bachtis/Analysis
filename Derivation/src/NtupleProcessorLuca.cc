#include "KaMuCa/Derivation/interface/NtupleProcessorLuca.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"

NtupleProcessorLuca::NtupleProcessorLuca(const std::string& outputFileName, bool isdata,bool gen) {

  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();

  isData=isdata;
  isGen=gen;
  

  w= new RooWorkspace("w","w");
  w->factory("c1[0,1]");
  w->factory("c2[0,1]");
  w->factory("gc1[0,1]");
  w->factory("gc2[0,1]");
  w->factory("eta1[-2.6,2.6]");
  w->factory("eta2[-2.6,2.6]");
  w->factory("phi1[-4,4]");
  w->factory("phi2[-4,4]");
  w->factory("cErr1[0,100]");
  w->factory("cErr2[0,100]");
  w->factory("mass[-200,200]");
  w->factory("massErr[0,200]");
  w->factory("N1[0,200]");
  w->factory("N2[0,200]");
  



  if (!(isdata||gen) )   
    w->defineSet("vars","c1,c2,gc1,gc2,eta1,eta2,phi1,phi2,cErr1,cErr2,mass,massErr,N1,N2");
  else
    w->defineSet("vars","c1,c2,eta1,eta2,phi1,phi2,cErr1,cErr2,mass,massErr,N1,N2");
  
  data=new RooDataSet("data","data",*w->set("vars"));

}



void NtupleProcessorLuca::processTree(const std::string&fileName,const std::string& cut) {
  TFile *fIn = new TFile(fileName.c_str());
  
  TTree *t = (TTree*)fIn->Get("ZTreeProducer");

  //reduce it!
  fOut->cd();

  TTree* reduced = t->CopyTree(cut.c_str());

  
  double pt1;
  double pt2;
  double ptErr1;
  double ptErr2;
  double mcpt1;
  double mcpt2;
  double eta1;
  double eta2;
  double phi1;
  double phi2;
  double mass;

  //  int N1=0;
  //int N2=0;

  Double_t covPlus[9];
  Double_t covMinus[9];


  reduced->SetBranchAddress("MuPos_pt",&pt1);
  reduced->SetBranchAddress("MuNeg_pt",&pt2);
  if (!isData) {
    reduced->SetBranchAddress("MuPosGenStatus1_pt",&mcpt1);
    reduced->SetBranchAddress("MuNegGenStatus1_pt",&mcpt2);
  }

  reduced->SetBranchAddress("MuPos_eta",&eta1);
  reduced->SetBranchAddress("MuNeg_eta",&eta2);
  reduced->SetBranchAddress("MuPos_phi",&phi1);
  reduced->SetBranchAddress("MuNeg_phi",&phi2);
  reduced->SetBranchAddress("Z_mass",&mass);
  reduced->SetBranchAddress("MuPosCovMatrix",covPlus);
  reduced->SetBranchAddress("MuNegCovMatrix",covMinus);



  for (int i=0;i<reduced->GetEntries();++i) {
    reduced->GetEntry(i);


    if (!isGen) {
      w->var("c1")->setVal(1.0/pt1);
      w->var("c2")->setVal(1.0/pt2);
      if (!isData) {
	w->var("gc1")->setVal(1.0/mcpt1);
	w->var("gc2")->setVal(1.0/mcpt2);
      }
      
    }
    else {
      w->var("c1")->setVal(1.0/mcpt1);
      w->var("c2")->setVal(1.0/mcpt2);
    }
    w->var("eta1")->setVal(eta1);
    w->var("eta2")->setVal(eta2);
    w->var("phi1")->setVal(phi1);
    w->var("phi2")->setVal(phi2);


    TLorentzVector a(1,1,1,1);
    a.SetPtEtaPhiM(pt1,eta1,phi1,0.10565);
    TLorentzVector b(1,1,1,1);
    b.SetPtEtaPhiM(pt2,eta2,phi2,0.10565);


    TMatrixDSym covMatrixPlus(3);
    covMatrixPlus(0,0) = covPlus[0];
    covMatrixPlus(0,1) = covPlus[1];
    covMatrixPlus(0,2) = covPlus[2];
    covMatrixPlus(1,0) = covPlus[3];
    covMatrixPlus(1,1) = covPlus[4];
    covMatrixPlus(1,2) = covPlus[5];
    covMatrixPlus(2,0) = covPlus[6];
    covMatrixPlus(2,1) = covPlus[7];
    covMatrixPlus(2,2) = covPlus[8];
    TMatrixD jacobianPlus(1,3);
    jacobianPlus(0,0)=a.Px()/a.Pt();
    jacobianPlus(0,1)=a.Py()/a.Pt();
    jacobianPlus(0,2)=0.0;
    TMatrixDSym dmCov = covMatrixPlus.Similarity(jacobianPlus);
    double error = dmCov(0,0) ;   
    if (error>=0)
      error=sqrt(error);
    else 
      error=0.0;
    ptErr1 = error/pt1;
    w->var("cErr1")->setVal(ptErr1);


    TMatrixDSym covMatrixMinus(3);
    covMatrixMinus(0,0) = covMinus[0];
    covMatrixMinus(0,1) = covMinus[1];
    covMatrixMinus(0,2) = covMinus[2];
    covMatrixMinus(1,0) = covMinus[3];
    covMatrixMinus(1,1) = covMinus[4];
    covMatrixMinus(1,2) = covMinus[5];
    covMatrixMinus(2,0) = covMinus[6];
    covMatrixMinus(2,1) = covMinus[7];
    covMatrixMinus(2,2) = covMinus[8];
    TMatrixD jacobianMinus(1,3);
    jacobianMinus(0,0)=b.Px()/b.Pt();
    jacobianMinus(0,1)=b.Py()/b.Pt();
    jacobianMinus(0,2)=0.0;
    dmCov = covMatrixMinus.Similarity(jacobianMinus);
    error = dmCov(0,0) ;   
    if (error>=0)
      error=sqrt(error);
    else 
      error=0.0;

    ptErr2 = error/pt2;
    w->var("cErr2")->setVal(ptErr2);


    w->var("N1")->setVal(0.0);
    w->var("N2")->setVal(0.0);

    if (!isGen) {
      w->var("mass")->setVal(mass);
    }
    else {
      TLorentzVector a(1,1,1,1);
      a.SetPtEtaPhiM(mcpt1,eta1,phi1,0.10565);
      TLorentzVector b(1,1,1,1);
      b.SetPtEtaPhiM(mcpt2,eta2,phi2,0.10565);
      w->var("mass")->setVal((a+b).M());
    }

    w->var("massErr")->setVal(0.5*sqrt(ptErr1*ptErr1+ptErr2*ptErr2)*mass);
    RooArgSet set(*w->set("vars"));
    data->add(set);
  }
  fIn->Close();

}



void NtupleProcessorLuca::write() {
  fOut->cd();
  data->Write();
  fOut->Close();
}
