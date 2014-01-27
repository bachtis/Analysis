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
            self.set('trackPt',shower.trackVector.Pt())
            self.set('trackEta',shower.trackVector.Eta())
            self.set('trackP',shower.trackMomentum)
            self.set('NECAL',float(len(shower.ecalConstituents)))
            self.set('NHCAL',float(len(shower.hcalConstituents)))
            self.fill()
            

