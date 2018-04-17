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
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

class MuonGenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
      explicit MuonGenAnalyzer(const edm::ParameterSet&);
      ~MuonGenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  

      // ----------member data ---------------------------
  edm::EDGetTokenT<reco::GenParticleCollection>      genParticles_;


  TFile *fout;
  TTree *tree;
  float mass;
  float rapidity;


};
MuonGenAnalyzer::MuonGenAnalyzer(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  genParticles_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));  
  fout=new TFile("muonTree.root","RECREATE");
  tree= new TTree("tree","tree");


  tree->Branch("mass",&mass,"mass/F");
  tree->Branch("rapidity",&rapidity,"rapidity/F");
}


MuonGenAnalyzer::~MuonGenAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}




//
// member functions
//

// ------------ method called for each event  ------------
void
MuonGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   reco::GenParticleCollection genMuons;

   Handle<reco::GenParticleCollection> genParticlesH;
   iEvent.getByToken(genParticles_,genParticlesH);

   for (unsigned int i=0;i<genParticlesH->size();++i) {
     const reco::GenParticle& p = (*genParticlesH)[i];
     if (abs(p.pdgId())==13 && p.status()==1 && fabs(p.eta())<2.4 && p.pt()>3.0)
       genMuons.push_back(p);
   }

   if (genMuons.size()<2)
     return;



   for (unsigned int i=0;i<genMuons.size()-1;++i) {
     for (unsigned int j=i+1;j<genMuons.size();++j) {
       reco::GenParticle mu1 = genMuons[i];
       reco::GenParticle mu2 = genMuons[j];
       if ((mu1.charge() + mu2.charge())==0 ) {
	 mass = (mu1.p4()+mu2.p4()).M();
	 rapidity = (mu1.p4()+mu2.p4()).Rapidity();
	 if (mass>3.0 && mass<3.2)
	   tree->Fill();
       }
     }
   }
}





// ------------ method called once each job just before starting event loop  ------------
void 
MuonGenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonGenAnalyzer::endJob() 
{
  fout->Write();
  fout->Close();

}

// ------------ method called when starting to processes a run  ------------
/*
void 
MuonGenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MuonGenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MuonGenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MuonGenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonGenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MuonGenAnalyzer);
