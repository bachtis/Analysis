from ROOT import TLorentzVector
import ROOT
from math import fabs,sqrt
from CMGTools.RootTools.utils.DeltaR import deltaR,deltaPhi
from array import array
import numpy


class ShowerFromChargedPion(object):
    def __init__(self,track,dr = 0.5):
        self.track=track
        self.dr = dr
        self.ecalEntrance =track.extrapolatedPoint(4).position()
        self.hcalEntrance =track.extrapolatedPoint(6).position()
        self.trackMomentum = self.track.trajectoryPoints()[0].momentum().energy()
        self.trackVector = self.track.trajectoryPoints()[0].momentum()
        self.ecalVector=None
        self.hcalVector=None
        self.ecalConstituents=[] 
        self.hcalConstituents=[] 
        self.hcalTiming=[]

        #for visualization
        #type = 0 track 1 ECAL , 2 HCAL

        self.typeO=[]
        self.eta=[]
        self.phi=[]
        self.rho=[]
        self.energy=[]


    def vectorFromConstituent(self,cluster):
        vec = ROOT.TVector3(cluster.position().x(),cluster.position().y(),cluster.position().z())
        vec= vec.Unit()
        vec*=cluster.energy()
        return ROOT.TLorentzVector(vec.x(),vec.y(),vec.z(),cluster.energy())

    def addConstituent(self,constituent,ecal = False):
        vec = self.vectorFromConstituent(constituent)

        if ecal and deltaR(vec.Eta(),vec.Phi(),self.trackVector.Eta(),self.trackVector.Phi())<self.dr:
            self.ecalConstituents.append(constituent)
            if self.ecalVector is None:
                self.ecalVector = vec
            else:    
                self.ecalVector += vec
            self.typeO.append(1)

            self.eta.append(constituent.position().Eta()-self.ecalEntrance.Eta())
            self.phi.append(deltaPhi(constituent.position().Phi(),self.ecalEntrance.Phi()))
            self.rho.append(constituent.position().Rho()-self.ecalEntrance.Rho())
            self.energy.append(constituent.energy())
            
            return True
        
        elif  deltaR(vec.Eta(),vec.Phi(),self.trackVector.Eta(),self.trackVector.Phi())<self.dr:
            self.hcalConstituents.append(constituent)
            self.hcalTiming.append(constituent.time())
            if self.hcalVector is None:
                self.hcalVector = vec
            else:    
                self.hcalVector += vec
            self.typeO.append(2)

            self.eta.append(constituent.position().Eta()-self.ecalEntrance.Eta())
            self.phi.append(deltaPhi(constituent.position().Phi(),self.ecalEntrance.Phi()))
            self.rho.append(constituent.position().Rho()-self.ecalEntrance.Rho())
            self.energy.append(constituent.energy())
            return True
        else:
            return False



    def calculate(self):
        self.hcalEtaPhi = ROOT.TH2F("hcalEtaPhi","",20,-0.5,0.5,-0.5,0.5)
        self.hcalEtaRho = ROOT.TH2F("hcalEtaRho","",20,-0.5,0.5,-0.5,0.5)

    def makeVisTree(self,tfile,treename):
        tfile.cd()
        tree =ROOT.TTree(treename,treename)

        typeO=numpy.zeros(1,float)
        tree.Branch("type",typeO,'type/D')

        eta=numpy.zeros(1,float)
        tree.Branch("eta",eta,'eta/D')

        phi=numpy.zeros(1,float)
        tree.Branch("phi",phi,'phi/D')

        rho=numpy.zeros(1,float)
        tree.Branch("rho",rho,'rho/D')

        energy=numpy.zeros(1,float)
        tree.Branch("energy",energy,'energy/D')

        for t,e,p,r,en in zip(self.typeO,self.eta,self.phi,self.rho,self.energy):
            typeO[0]=t
            eta[0]=e
            phi[0]=p
            rho[0]=r
            energy[0]=en

            tree.Fill()
        tree.Write()
        
        

    
    def __str__(self):
        str = """
        New Shower from Pion
        --------------------
        """
        str=str+ 'track pt,eta,phi,p,ecal eta, ecal phi, hcal eta, hcal phi  {pt},{eta},{phi},{p},{eeta},{ephi},{heta},{hphi} '.format(pt=self.trackVector.Pt(),eta=self.trackVector.Eta(),phi=self.trackVector.Phi(),p=self.trackMomentum,eeta=self.ecalEntrance.Eta(),ephi=self.ecalEntrance.Phi(),heta=self.hcalEntrance.Eta(),hphi=self.hcalEntrance.Phi())+'\n'
        str=str+  """ECAL constituents
                     -----------------
                     """
        for c in self.ecalConstituents:
            str=str+ 'ID={id},Drho={rho},eta={eta},Dphi={phi},z={z},energy={energy} \n'.format(id= c.detId(),rho=c.position().rho(),eta=c.position().eta(),phi=c.position().phi(),z=c.position().z(),energy=c.energy())

        str=str+  """HCAL constituents
                     -----------------
                     """
        for c in self.hcalConstituents:
            str=str+ 'ID={id},Drho={rho},Deta={eta},Dphi={phi},z={z},energy={energy}\n'.format(id= c.detId(),rho=c.position().rho(),eta=c.position().eta(),phi=c.position().phi(),z=c.position().z(),energy=c.energy())


        return str
            
