// -*- C++ -*-
//
// Package:    SpikeFix
// Class:      SpikeFix
// 
/**\class SpikeFix SpikeFix.cc PFKludge/SpikeFix/plugins/SpikeFix.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Michail Bachtis
//         Created:  Wed, 10 Sep 2014 15:45:23 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Math/interface/deltaR.h"
//
// class declaration
//

class SpikeFix : public edm::EDProducer {
   public:
      explicit SpikeFix(const edm::ParameterSet&);
      ~SpikeFix();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SpikeFix::SpikeFix(const edm::ParameterSet& iConfig)
{
  produces<reco::PFCandidateCollection>();
  
}


SpikeFix::~SpikeFix()
{

}

void
SpikeFix::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<reco::PFCandidateCollection> pfH;
   iEvent.getByLabel("particleFlow",pfH);

   Handle<reco::PFJetCollection> pfJH;
   iEvent.getByLabel("ak5PFJets",pfJH);

   Handle<reco::CaloJetCollection> caloJH;
   iEvent.getByLabel("ak5CaloJets",caloJH);

   std::auto_ptr<reco::PFCandidateCollection> out(new reco::PFCandidateCollection);


   for (auto cand: *pfH) {
     
     if (cand.particleId()==5 && fabs(cand.eta()) >2.85  && fabs(cand.eta())<2.95) {
       //find the nearest PF jet and calojet
       float hcalPF=0.0;
       float hfPF=0.0;
       float hcalCalo=0.0;
       float hfCalo=0.0;
       
       double dr = 1000.0;

       for (auto jet : *pfJH) {
	 double delta = deltaR(cand.eta(),cand.phi(),jet.eta(),jet.phi());
	 if (delta<0.1 && delta<dr) {
	   dr=delta;
	   hcalPF = jet.neutralHadronEnergy();
	   hfPF = jet.HFHadronEnergy();
	 }
       }
       dr=1000.0;
       for (auto jet : *caloJH) {
	 double delta = deltaR(cand.eta(),cand.phi(),jet.eta(),jet.phi());
	 if (delta<0.1 && delta<dr) {
	   dr=delta;
	   hcalCalo = jet.hadEnergyInHE();
	   hfCalo = jet.hadEnergyInHF();
	 }
       }

       float scaleFactor =1.0;

       if(hcalPF>0.0) {
	 scaleFactor = (hcalCalo+hfCalo-hfPF)/hcalPF;
       }

       reco::PFCandidate c = cand;
       c.rescaleMomentum(scaleFactor);
       out->push_back(c);
     }
     else {
       out->push_back(cand);
     }

   }

   iEvent.put(out);

}

// ------------ method called once each job just before starting event loop  ------------
void 
SpikeFix::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SpikeFix::endJob() {
}

void
SpikeFix::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SpikeFix);
