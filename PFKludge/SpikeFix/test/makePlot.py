import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle



jetH = Handle('std::vector<reco::PFJet>')
jetH2 = Handle('std::vector<reco::PFJet>')
jetH3 = Handle('std::vector<reco::CaloJet>')
pfH = Handle('std::vector<reco::PFCandidate>')
pfH2 = Handle('std::vector<reco::PFCandidate>')





events = Events([
'/tmp/bachtis/out.root'
]
)


etaJet    = ROOT.TH1F("etaJet","",100,-5.,5.)
etaJet2    = ROOT.TH1F("etaJet2","",100,-5.,5.)
etaJet3    = ROOT.TH1F("etaJet3","",100,-5.,5.)

eta    = ROOT.TH1F("eta","",100,-5.,5.)
eta2    = ROOT.TH1F("eta2","",100,-5.,5.)


for event in events:
    event.getByLabel('ak5PFJets','','RECO',jetH)
    event.getByLabel('ak5PFJets','','FIX',jetH2)
    event.getByLabel('ak5CaloJets',jetH3)
    event.getByLabel('particleFlow',pfH)
    event.getByLabel('fixedPF',pfH2)

    jet = jetH.product()
    jet2 = jetH2.product()
    jet3 = jetH3.product()
    pf = pfH.product()
    pf2 = pfH2.product()
    
    for j in jet:
        if j.pt()>30:
            etaJet.Fill(j.eta())
    for j in jet2:
        if j.pt()>30:
            etaJet2.Fill(j.eta())

    for j in jet3:
        if j.pt()>30:
            etaJet3.Fill(j.eta())


    for p in pf:
        if p.charge()==0 and p.pdgId() !=22:
            eta.Fill(p.eta(),p.energy())

    for p in pf2:
        if p.charge()==0 and p.pdgId() !=22:
            eta2.Fill(p.eta(),p.energy())



etaJet.Draw();
etaJet2.SetLineColor(ROOT.kRed)
etaJet2.Draw("SAME");

l=ROOT.TLegend(0.4,0.6,0.8,0.9)
l.AddEntry(etaJet,"before","l")
l.AddEntry(etaJet2,"after","l")
l.Draw()
