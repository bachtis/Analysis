// -*- C++ -*-
//
// Package:    TrackAnalysis/HitAnalyzer
// Class:      HitAnalyzer
// 
/**\class HitAnalyzer HitAnalyzer.cc TrackAnalysis/HitAnalyzer/plugins/HitAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michail Bachtis
//         Created:  Mon, 21 Mar 2016 14:17:37 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

class MuonCalibAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonCalibAnalyzer(const edm::ParameterSet&);
      ~MuonCalibAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  bool muonIDZ(const pat::Muon&,const reco::Vertex&);
  bool muonIDOnia(const pat::Muon&,const reco::Vertex&);
  bool selectResonance(const pat::Muon&, const pat::Muon&,const reco::Vertex&);
  int  findGenParticle(const pat::Muon&, const pat::PackedGenParticleCollection&);
  

      // ----------member data ---------------------------
  edm::EDGetTokenT<pat::MuonCollection>      muons_;
  edm::EDGetTokenT<pat::PackedGenParticleCollection>      genParticles_;
  edm::EDGetTokenT<reco::VertexCollection>      vertices_;
  edm::EDGetTokenT<pat::METCollection>      mets_;


  TFile *fout;
  TTree *tree;

  float pt1;
  float c1;
  float sc1;
  float dsc1;

  float gc1;
  float eta1;
  float phi1;
  float dpt1;
  float dxy1;
  float dz1;
  float mcpt1;

  float pt2;
  float c2;
  float sc2;
  float dsc2;
  float gc2;
  float eta2;
  float phi2;
  float dpt2;
  float dxy2;
  float dz2;
  float mcpt2;

  float mass;
  float massErr;
  float genMass;
  float recoil;

  int run;

  bool isOnia_;

};

bool MuonCalibAnalyzer::muonIDZ(const pat::Muon& muon, const reco::Vertex& vertex) {
  bool acc = muon.pt()>10&&fabs(muon.eta())<2.4;
  bool id = muon::isTightMuon(muon,vertex);
  bool iso = muon.isolationR03().sumPt/muon.pt()<0.2;
  bool vtx = fabs(muon.innerTrack()->dxy(vertex.position()))<0.02 &&  fabs(muon.innerTrack()->dz(vertex.position()))<0.2;

  return id&&vtx&&iso&&acc;

}


bool MuonCalibAnalyzer::muonIDOnia(const pat::Muon& muon, const reco::Vertex& vertex) {
  bool acc = muon.pt()>3.0&&fabs(muon.eta())<2.4;
  bool id = muon::isSoftMuon(muon,vertex);
  bool vtx = fabs(muon.innerTrack()->dxy(vertex.position()))<0.02 &&  fabs(muon.innerTrack()->dz(vertex.position()))<0.2;
  return id&&vtx&&acc;
}

bool MuonCalibAnalyzer::selectResonance(const pat::Muon& m1, const pat::Muon& m2,const reco::Vertex& vertex) {

  if (!(m1.isTrackerMuon()&& m2.isTrackerMuon()))
    return false;



  if (!(m1.charge()+m2.charge()==0))
    return false;


  if (isOnia_) {
    if (!muonIDOnia(m1,vertex))
      return false;
    if (!muonIDOnia(m2,vertex))
      return false;

    float m = (m1.p4()+m2.p4()).M();

    if ((m>2.7&&m<3.5) || (m>9&&m<10))
      return true;

    return false;

  }
  else {
    if (!muonIDZ(m1,vertex))
      return false;
    if (!muonIDZ(m2,vertex))
      return false;

    float m = (m1.p4()+m2.p4()).M();

    if (m>70&&m<120)
      return true;

    return false;

  }


}


int  MuonCalibAnalyzer::findGenParticle(const pat::Muon& mu,const pat::PackedGenParticleCollection& gen) {
  float drMin=10000;
   int index=-1;
  for (unsigned int i=0;i<gen.size();++i) {
    const pat::PackedGenParticle& p=gen[i];
    float dr=deltaR(p.eta(),p.phi(),mu.eta(),mu.phi());
    if (dr<0.1) {
      if (dr<drMin) {
	drMin=dr;
	index=(int)i;
      }

    }

  }  
  return index;
}





//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonCalibAnalyzer::MuonCalibAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  muons_ = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"));
  genParticles_ = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));  
  vertices_ = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));  
  mets_ = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"));  

  isOnia_ = iConfig.getParameter<bool>("isOnia");

  fout=new TFile("muonTree.root","RECREATE");
  tree= new TTree("tree","tree");


  tree->Branch("pt1",&pt1,"pt1/F");
  tree->Branch("pt2",&pt2,"pt2/F");
  tree->Branch("mcpt1",&mcpt1,"mcpt1/F");
  tree->Branch("mcpt2",&mcpt2,"mcpt2/F");
  tree->Branch("gc1",&gc1,"gc1/F");
  tree->Branch("gc2",&gc2,"gc2/F");
  tree->Branch("c1",&c1,"c1/F");
  tree->Branch("c2",&c2,"c2/F");
  tree->Branch("sc1",&sc1,"sc1/F");
  tree->Branch("sc2",&sc2,"sc2/F");
  tree->Branch("dsc1",&dsc1,"dsc1/F");
  tree->Branch("dsc2",&dsc2,"dsc2/F");

  tree->Branch("eta1",&eta1,"eta1/F");
  tree->Branch("eta2",&eta2,"eta2/F");
  tree->Branch("phi1",&phi1,"phi1/F");
  tree->Branch("phi2",&phi2,"phi2/F");
  tree->Branch("cErr1",&dpt1,"cErr1/F");
  tree->Branch("cErr2",&dpt2,"cErr2/F");
  tree->Branch("dxy1",&dxy1,"dxy1/F");
  tree->Branch("dxy2",&dxy2,"dxy2/F");
  tree->Branch("dz1",&dz1,"dz1/F");
  tree->Branch("dz2",&dz2,"dz2/F");
  
  tree->Branch("mass",&mass,"mass/F");
  tree->Branch("massErr",&massErr,"massErr/F");
  tree->Branch("genMass",&genMass,"genMass/F");
  tree->Branch("run",&run,"run/I");
  tree->Branch("recoil",&recoil,"recoil/F");
  
}


MuonCalibAnalyzer::~MuonCalibAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}




//
// member functions
//

// ------------ method called for each event  ------------
void
MuonCalibAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   pat::PackedGenParticleCollection genMuons;


   if (iEvent.run()<50000) {
     Handle<pat::PackedGenParticleCollection> genParticlesH;
     iEvent.getByToken(genParticles_,genParticlesH);
     for (unsigned int i=0;i<genParticlesH->size();++i) {
       const pat::PackedGenParticle& p = (*genParticlesH)[i];
       if (abs(p.pdgId())==13 and p.status()==1)
	 genMuons.push_back(p);
     }
     
   }

   run=iEvent.run();

   Handle<reco::VertexCollection> vertexH;
   iEvent.getByToken(vertices_,vertexH);
   
   if (vertexH->size()==0)
     return;
       

   const reco::Vertex& vertex = vertexH->at(0);
   

   Handle<pat::MuonCollection> muonH;
   iEvent.getByToken(muons_,muonH);


   unsigned N = muonH->size();

   if (N<2)
     return;



   Handle<pat::METCollection> metH;
   iEvent.getByToken(mets_,metH);

   reco::Candidate::LorentzVector metV = (*metH)[0].p4();
   

   for (unsigned int i=0;i<N-1;++i) {
     for (unsigned int j=i+1;j<N;++j) {

       const pat::Muon& mu1  = (*muonH)[i]; 
       const pat::Muon& mu2  = (*muonH)[j];        


       if (selectResonance((*muonH)[i],(*muonH)[j],vertex)) {
	   const pat::Muon& pos  = mu1.charge()>0 ? mu1 : mu2;
	   const pat::Muon& neg  = mu1.charge()<0 ? mu1 : mu2;


	   recoil=(-(pos.p4()+neg.p4()+metV)).pt();
   
	   pt1=pos.innerTrack()->pt();
	   if (pos.isGlobalMuon()&&pos.standAloneMuon().isNonnull()) {
	     sc1=1.0/pos.standAloneMuon()->pt();
	     dsc1=pos.standAloneMuon()->ptError()/pos.standAloneMuon()->pt();
	   }
	   else {
	     sc1=-1.0;
	     dsc1=-1.0;
	   }


	   dpt1=pos.innerTrack()->ptError()/pt1;
	   c1=1.0/pt1;
	   eta1=pos.innerTrack()->eta();
	   phi1=pos.innerTrack()->phi();
	   int ptr1 = findGenParticle(pos,genMuons);

	   //	   printf("Track length=%f\n",pos.innerTrack()->outerPosition().rho()-pos.innerTrack()->innerPosition().rho());
	   
	   if (ptr1>=0)  {
	     mcpt1=genMuons[ptr1].pt();
	     gc1=1.0/mcpt1;
	   }
	   else {
	     mcpt1=-1.0;
	     gc1=-1.0;
	   }

	   dxy1=pos.innerTrack()->dxy(vertex.position());
	   dxy2=neg.innerTrack()->dxy(vertex.position());
	   dz1=pos.innerTrack()->dz(vertex.position());
	   dz2=neg.innerTrack()->dz(vertex.position());


	   pt2=neg.innerTrack()->pt();

	   if (neg.isGlobalMuon()&&neg.standAloneMuon().isNonnull()) {
	     sc2=1.0/neg.standAloneMuon()->pt();
	     dsc2=neg.standAloneMuon()->ptError()/neg.standAloneMuon()->pt();
	   }


	   dpt2=neg.innerTrack()->ptError()/pt2;
	   c2=1.0/pt2;
	   eta2=neg.innerTrack()->eta();
	   phi2=neg.innerTrack()->phi();
	   int ptr2 = findGenParticle(neg,genMuons);


	   if (ptr2>=0)  {
	     mcpt2=genMuons[ptr2].pt();
	     gc2=1.0/mcpt2;
	   }
	   else {
	     mcpt2=-1;
	     gc2=-1;
	   }
	   
	   math::PtEtaPhiMLorentzVector posVec(pt1,eta1,phi1,0.105658);
	   math::PtEtaPhiMLorentzVector negVec(pt2,eta2,phi2,0.105658);
	   mass=(posVec+negVec).M();
	   massErr = 0.5*mass*sqrt(dpt1*dpt1+dpt2*dpt2);

	   genMass=-1.0;
	   if (iEvent.run()<50000 && mcpt1>0 && mcpt2>0) {
	     math::PtEtaPhiMLorentzVector posVec(mcpt1,eta1,phi1,0.105658);
	     math::PtEtaPhiMLorentzVector negVec(mcpt2,eta2,phi2,0.105658);
	     genMass = (posVec+negVec).M();
	   }

	   tree->Fill();
	   return;
	 }
     }
   }

}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonCalibAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonCalibAnalyzer::endJob() 
{
  fout->Write();
  fout->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonCalibAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonCalibAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonCalibAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonCalibAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonCalibAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonCalibAnalyzer);
