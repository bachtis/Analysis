#include "KaMuCa/Derivation/interface/NtupleProcessor.h"
#include "TLorentzVector.h"

NtupleProcessor::NtupleProcessor(const std::string& outputFileName, bool isdata,bool gen,const std::string& treePrefix) {

  fOut = new TFile(outputFileName.c_str(),"RECREATE");
  fOut->cd();

  isData=isdata;
  isGen=gen;
  prefix=treePrefix;

  

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



void NtupleProcessor::processTree(const std::string&fileName,const std::string& cut) {
  TFile *fIn = new TFile(fileName.c_str());


  
  TTree *t = (TTree*)fIn->Get("tree");

  //reduce it!
  fOut->cd();
  TTree* reduced = t->CopyTree(cut.c_str());
  
  float pt1;
  float pt2;
  float ptErr1;
  float ptErr2;
  float mcpt1;
  float mcpt2;
  float eta1;
  float eta2;
  float phi1;
  float phi2;
  float mass;

  int N1=0;
  int N2=0;
  int q1=0;
  int q2=0;


  reduced->SetBranchAddress((prefix+"_l1_charge").c_str(),&q1);
  reduced->SetBranchAddress((prefix+"_l2_charge").c_str(),&q2);

  reduced->SetBranchAddress((prefix+"_l1_pt").c_str(),&pt1);
  reduced->SetBranchAddress((prefix+"_l2_pt").c_str(),&pt2);
  if (!isData) {
    reduced->SetBranchAddress((prefix+"_l1_mcPt").c_str(),&mcpt1);
    reduced->SetBranchAddress((prefix+"_l2_mcPt").c_str(),&mcpt2);
  }
  reduced->SetBranchAddress((prefix+"_l1_ptErr").c_str(),&ptErr1);
  reduced->SetBranchAddress((prefix+"_l2_ptErr").c_str(),&ptErr2);
  reduced->SetBranchAddress((prefix+"_l1_eta").c_str(),&eta1);
  reduced->SetBranchAddress((prefix+"_l2_eta").c_str(),&eta2);
  reduced->SetBranchAddress((prefix+"_l1_phi").c_str(),&phi1);
  reduced->SetBranchAddress((prefix+"_l2_phi").c_str(),&phi2);
  reduced->SetBranchAddress((prefix+"_mll").c_str(),&mass);
  reduced->SetBranchAddress((prefix+"_l1_trackerHits").c_str(),&N1);
  reduced->SetBranchAddress((prefix+"_l2_trackerHits").c_str(),&N2);


  for (int i=0;i<reduced->GetEntries();++i) {
    reduced->GetEntry(i);


    if (q1>0) {

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
      w->var("cErr1")->setVal(ptErr1/pt1);
      w->var("cErr2")->setVal(ptErr2/pt2);
      w->var("N1")->setVal(720.0/(N1+4));
      w->var("N2")->setVal(720.0/(N2+4));

    }

    else {

      if (!isGen) {
	w->var("c2")->setVal(1.0/pt1);
	w->var("c1")->setVal(1.0/pt2);
	if (!isData) {
	  w->var("gc2")->setVal(1.0/mcpt1);
	  w->var("gc1")->setVal(1.0/mcpt2);
	}

      }
      else {
	w->var("c2")->setVal(1.0/mcpt1);
	w->var("c1")->setVal(1.0/mcpt2);
      }
      w->var("eta2")->setVal(eta1);
      w->var("eta1")->setVal(eta2);
      w->var("phi2")->setVal(phi1);
      w->var("phi1")->setVal(phi2);
      w->var("cErr2")->setVal(ptErr1/pt1);
      w->var("cErr1")->setVal(ptErr2/pt2);
      w->var("N2")->setVal(720.0/(N1+4));
      w->var("N1")->setVal(720.0/(N2+4));
    }
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



    w->var("massErr")->setVal(0.5*sqrt(ptErr1*ptErr1/(pt1*pt1)+ptErr2*ptErr2/(pt2*pt2))*mass);
    RooArgSet set(*w->set("vars"));
    data->add(set);

  }

  fIn->Close();

}



void NtupleProcessor::write() {
  fOut->cd();
  data->Write();
  fOut->Close();
}
