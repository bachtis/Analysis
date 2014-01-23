from ROOT import TLorentzVector
import ROOT
from math import fabs,sqrt
from CMGTools.RootTools.utils.DeltaR import deltaR,deltaPhi

class ShowerFromChargedPion(object):
    def __init__(self,track,dr = 0.5):
        self.track=track
        self.dr = dr
        self.ecalEntrance =track.extrapolatedPoint(4).position()
        self.hcalEntrance =track.extrapolatedPoint(5).position()
        self.trackMomentum = self.track.extrapolatedPoint(4).momentum().energy()
        self.trackVector = self.track.extrapolatedPoint(1).momentum()
        self.ecalVector=None
        self.hcalVector=None
        self.ecalConstituents=[] 
        self.hcalConstituents=[] 
        self.hcalTiming=[]


    def vectorFromConstituent(self,cluster):
        vec = ROOT.TVector3(cluster.position().x(),cluster.position().y(),cluster.position().z())
        vec= vec.Unit()
        vec*=cluster.energy()
        return ROOT.TLorentzVector(vec.x(),vec.y(),vec.z(),cluster.energy())

    def addConstituent(self,constituent,ecal = False):
        vec = self.vectorFromConstituent(constituent)

        if ecal and deltaR(vec.Eta(),vec.Phi(),self.ecalEntrance.Eta(),self.ecalEntrance.Phi())<self.dr:
            self.ecalConstituents.append(constituent)
            if self.ecalVector is None:
                self.ecalVector = vec
            else:    
                self.ecalVector += vec
            return True
        
        elif  deltaR(vec.Eta(),vec.Phi(),self.ecalEntrance.Eta(),self.ecalEntrance.Phi())<self.dr:
            self.hcalConstituents.append(constituent)
            self.hcalTiming.append(constituent.time())
            if self.hcalVector is None:
                self.hcalVector = vec
            else:    
                self.hcalVector += vec

            return True
        else:
            return False



    def calculate(self):
        self.hcalEtaPhi = ROOT.TH2F("hcalEtaPhi","",20,-0.5,0.5,-0.5,0.5)
        self.hcalEtaRho = ROOT.TH2F("hcalEtaRho","",20,-0.5,0.5,-0.5,0.5)
        

    def __str__(self):
        str = """
        New Shower from Pion
        --------------------
        """
        str=str+ 'track pt,eta,phi,p {pt},{eta},{phi},{p} '.format(pt=self.trackVector.Pt(),eta=self.trackVector.Eta(),phi=self.trackVector.Phi(),p=self.trackMomentum)+'\n'
        str=str+  """ECAL constituents
                     -----------------
                     """
        for c in self.ecalConstituents:
            str=str+ 'ID={id},Drho={rho},Deta={eta},Dphi={phi},z={z},energy={energy}\n'.format(id= c.detId(),rho=c.position().rho()-self.ecalEntrance.Rho(),eta=c.position().eta()-self.ecalEntrance.Eta(),phi=c.position().phi()-self.ecalEntrance.Phi(),z=c.position().z(),energy=c.energy())

        str=str+  """HCAL constituents
                     -----------------
                     """
        for c in self.hcalConstituents:
            str=str+ 'ID={id},Drho={rho},Deta={eta},Dphi={phi},z={z},energy={energy}\n'.format(id= c.detId(),rho=c.position().rho()-self.hcalEntrance.Rho(),eta=c.position().eta()-self.hcalEntrance.Eta(),phi=c.position().phi()-self.hcalEntrance.Phi(),z=c.position().z(),energy=c.energy())


        return str
            
