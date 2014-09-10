import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle



jetH = Handle('std::vector<reco::PFJet>')
jetH2 = Handle('std::vector<reco::PFJet>')





events = Events([
'/tmp/bachtis/out.root'
]
)


etaJet    = ROOT.TH1F("etaJet","",100,-5.,5.)
etaJet2    = ROOT.TH1F("etaJet2","",100,-5.,5.)


for event in events:
    event.getByLabel('ak5PFJets','','RECO',jetH)
    event.getByLabel('ak5PFJets','','FIX',jetH2)

    jet = jetH.product()
    jet2 = jetH2.product()

    
    for j in jet:
        etaJet.Fill(j.eta())
    for j in jet2:
        etaJet2.Fill(j.eta())




etaJet.Draw();
etaJet2.SetLineColor(ROOT.kRed)
etaJet2.Draw("SAME");

l=ROOT.TLegend(0.4,0.6,0.8,0.9)
l.AddEntry(etaJet,"before","l")
l.AddEntry(etaJet2,"after","l")
l.Draw()
