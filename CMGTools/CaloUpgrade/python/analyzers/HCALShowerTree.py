import operator 
import itertools
import copy
from math import fabs
from ROOT import TLorentzVector

from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy
        
class HCALShowerTree( TreeAnalyzerNumpy ):

    def declareVariables(self):

        self.branch('nPU')
        self.branch('trackPt')
        self.branch('trackP')
        self.branch('trackEta')
        self.branch('NECAL')
        self.branch('NHCAL')
        
        super(HCALShowerTree, self).declareVariables()
        
 
       
    def process(self, iEvent, event):
        for shower in event.showers:
            self.reset()
            self.set('nPU',event.puInteractions)
            self.set('trackPt',block.trackVector.Pt())
            self.set('trackEta',block.trackVector.Eta())
            self.set('trackP',block.trackMomentum)
            self.set('NECAL',float(len(block.ecalConstituents)))
            self.set('NHCAL',float(len(block.hcalConstiturents)))
            self.fill()
            

