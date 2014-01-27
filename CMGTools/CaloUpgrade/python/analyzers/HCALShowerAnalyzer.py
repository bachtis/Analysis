import operator 
import itertools
import copy
from math import fabs,sqrt
from sets import Set

from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.Event import Event
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.CaloUpgrade.tools.DataFormats import ShowerFromChargedPion 
import ROOT

        
        
class HCALShowerAnalyzer( Analyzer ):

    def __init__(self, cfg_ana, cfg_comp, looperName ):
        self.doVis=True
        self.visFile=ROOT.TFile("visInput.root","RECREATE")
        self.doVis=True
        self.counter=0
        super(HCALShowerAnalyzer,self).__init__(cfg_ana,cfg_comp,looperName)



        
    def declareHandles(self):
        ''' Here declare handles of all objects we possibly need
        '''
        super(HCALShowerAnalyzer, self).declareHandles()

        self.handles['hcalHits'] = AutoHandle( ('particleFlowRecHitHCAL','',''),'std::vector<reco::PFRecHit>')
        self.handles['hoHits'] = AutoHandle( ('particleFlowRecHitHO',''),'std::vector<reco::PFRecHit>')
        self.handles['ecalHits'] = AutoHandle( ('particleFlowRecHitECAL',''),'std::vector<reco::PFRecHit>')
        self.handles['puInfo'] = AutoHandle( ('addPileupInfo',''),'std::vector<PileupSummaryInfo>')
        self.handles['tracks'] = AutoHandle( ('pfTrack',''),'std::vector<reco::PFRecTrack>')

    def beginLoop(self):
        super(HCALShowerAnalyzer,self).beginLoop()
        
       
    def process(self, iEvent, event):
        self.event = iEvent.eventAuxiliary().id().event()
        self.readCollections( iEvent )
        event.pu = self.handles['puInfo'].product()
        event.puInteractions=event.pu[0].getTrueNumInteractions()

        event.tracks = self.handles['tracks'].product()
        event.hcalHits = self.handles['hcalHits'].product()
        event.ecalHits = self.handles['ecalHits'].product()
        event.hoHits = self.handles['hoHits'].product()

        event.showers=[]

        #loop on tracks and match them with clusters
        for track in event.tracks:
            shower = ShowerFromChargedPion(track,0.5)
            for hit in event.ecalHits:
                shower.addConstituent(hit,True)
            for hit in (event.hcalHits):
                shower.addConstituent(hit)
            for hit in (event.hoHits):
                shower.addConstituent(hit)
            event.showers.append(shower)    
            print shower


        if self.doVis:
            for i,shower in enumerate(event.showers):
                shower.makeVisTree(self.visFile,"t_"+str(self.counter)+"_"+str(i))

        self.counter=self.counter+1
        return True
    

