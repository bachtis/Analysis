import ROOT
import itertools
import math
from DataFormats.FWLite import Events, Handle

ROOT.gROOT.ProcessLine(".x tdrstyle.C")


def deltaPhi( p1, p2):
    '''Computes delta phi, handling periodic limit conditions.'''
    res = p1 - p2
    while res > math.pi:
        res -= 2*math.pi
    while res < -math.pi:
        res += 2*math.pi
    return res


events = Events ('patTuple_standard.root')
muonHandle  = Handle ('std::vector<pat::Muon>')
vertexHandle  = Handle ('std::vector<reco::Vertex>')
pfWHandle  = Handle ('std::vector<reco::PFCandidate>')
pfHandle  = Handle ('std::vector<reco::PFCandidate>')

def add(collection):
    sum=0
    sumV=collection[0].p4()
    for i,object in enumerate(collection):
            sum=sum+object.pt()
            if i==0:
                continue
            sumV=sumV+object.p4()
    return sum,sumV;        

profile = ROOT.TProfile("sumET","sumET",25,0,50,0,1000.)
profile.GetXaxis().SetTitle("# vertices")
profile.GetYaxis().SetTitle("<neutral #Sigma E_{T}>")
profile.SetMarkerColor(ROOT.kBlue)

profileW = ROOT.TProfile("sumETW","sumETW",25,0,50,0,1000.)
profileW.GetXaxis().SetTitle("# vertices")
profileW.GetYaxis().SetTitle("<neutral #Sigma E_{T}>")
profileW.SetMarkerColor(ROOT.kBlue)

histoMET = ROOT.TH1F("pfMET","sumET",20,0,100)
histoMETW = ROOT.TH1F("pfMETW","sumET",20,0,100)





## weightV = ROOT.TProfile("weightV","weight",25,0,50,0,1.)
## weightV.GetXaxis().SetTitle("# vertices")
## weightV.GetYaxis().SetTitle("<w_{i}>")

## weightEta = ROOT.TProfile("weightEta","weight",25,0,5,0,1.)
## weightEta.GetXaxis().SetTitle("#eta")
## weightEta.GetYaxis().SetTitle("<w_{i}>")

## weightPhi = ROOT.TProfile("weightPhi","weight",32,-3.14,3.14,0,1.)
## weightPhi.GetXaxis().SetTitle("#phi")
## weightPhi.GetYaxis().SetTitle("<w_{i}>")


## weightE = ROOT.TProfile("weightE","weight",25,0,100,0,1.)
## weightE.GetXaxis().SetTitle("E [GeV]")
## weightE.GetYaxis().SetTitle("<w_{i}>")




for event in events:
    event.getByLabel('selectedPatMuons',muonHandle)
    event.getByLabel('offlinePrimaryVertices',vertexHandle)
    event.getByLabel('particleFlowWeighted',pfWHandle)
    event.getByLabel('particleFlow',pfHandle)

    muons =muonHandle.product()
    vertex =vertexHandle.product()[0]
    vertices = len(vertexHandle.product())
    pf = pfHandle.product()
    pfW = pfWHandle.product()
    
##     for n in neutrals:
##         weightV.Fill(vertices,n.mva_e_pi())
##         weightEta.Fill(abs(n.eta()),n.mva_e_pi())
##         weightPhi.Fill(n.phi(),n.mva_e_pi())
##         weightE.Fill(n.energy(),n.mva_e_pi())

    sum,sumV=add(pf)
    sumW,sumVW=add(pfW)
    
    
    profile.Fill(vertices,sum)
    profileW.Fill(vertices,sumW)
    histoMET.Fill(sumV.rho())
    histoMETW.Fill(sumVW.rho())
    
    
