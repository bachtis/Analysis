import ROOT
import numpy
import math

class MCClosureTools (object):


    def loop(self,sign,data):
        h = ROOT.TH1F("res","res",500,-0.005,0.005)
        
        for evt in range(0,data.numEntries()):
            line = data.get(evt)
            if sign == 'pos':
                curv = 'curvRaw1'
                gencurv = 'curvGenRaw1'
            else:    
                curv = 'curvRaw2'
                gencurv = 'curvGenRaw2'
                
            h.Fill((1./line.find(gencurv).getVal()-1./line.find(curv).getVal())*line.find(gencurv).getVal())
        return h

    def loopMass(self,sign,data):
        h = ROOT.TH1F("res","res",500,0.97,1.03)
        
        for evt in range(0,data.numEntries()):
            line=data.get(evt)
            if line.find('curvGenRaw1').getVal()==0 or line.find('curvGenRaw2').getVal()==0:
                continue

            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvGenRaw1').getVal(),
                                   line.find('etaRaw1').getVal(),
                                   line.find('phiRaw1').getVal(),
                                   0.1056583715)
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvGenRaw2').getVal(),
                                   line.find('etaRaw2').getVal(),
                                   line.find('phiRaw2').getVal(),
                                   0.1056583715)
            h.Fill(line.find("massRaw").getVal()/(v1+v2).M())
        return h

    def loopMassProfile(self,sign,data):
        h = ROOT.TProfile("res","res",10,80,100,-5,5,'s')
        
        for evt in range(0,data.numEntries()):
            line=data.get(evt)
            if line.find('curvGenRaw1').getVal()==0 or line.find('curvGenRaw2').getVal()==0:
                continue

            v1=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(1./line.find('curvGenRaw1').getVal(),
                                   line.find('etaRaw1').getVal(),
                                   line.find('phiRaw1').getVal(),
                                   0.1056583715)
            
            v2=ROOT.TLorentzVector()
            v2.SetPtEtaPhiM(1./line.find('curvGenRaw2').getVal(),
                                   line.find('etaRaw2').getVal(),
                                   line.find('phiRaw2').getVal(),
                                   0.1056583715)
            h.Fill((v1+v2).M(),line.find("massRaw").getVal()-(v1+v2).M())
        return h
    

    def getMeanAndSpreadArithmetic(self,sign,data):
        h = self.loop(sign,data)
        hist=h
        return h.GetMean(),h.GetMeanError(),h.GetRMS(),h.GetRMSError()
    def getMeanAndSpreadArithmeticMass(self,sign,data):
        h = self.loopMass(sign,data)
        hist=h
        return h.GetMean(),h.GetMeanError(),h.GetRMS(),h.GetRMSError()

    def getMeanAndSpreadFit(self,sign,data):
        h= self.loop(sign,data)

        w=ROOT.RooWorkspace("w","w")
        w.factory("res[-0.003,0.003]")
        datahist = ROOT.RooDataHist("histD","histD",ROOT.RooArgList(w.var("res")),h)
        w.factory("RooGaussian::gaus(res,mean[-0.003,0.003],sigma[0,0.01])")
        w.pdf("gaus").fitTo(datahist)
        
        return w.var('mean').getVal(),w.var('mean').getError(),w.var('sigma').getVal(),w.var('sigma').getError()

    def getMeanAndSpreadFitMass(self,sign,data):
        h= self.loopMass(sign,data)

        w=ROOT.RooWorkspace("w","w")
        w.factory("res[0.97,1.03]")
        datahist = ROOT.RooDataHist("histD","histD",ROOT.RooArgList(w.var("res")),h)
        w.factory("RooGaussian::gaus(res,mean[0.97,1.03],sigma[0.,5.])")
        w.pdf("gaus").fitTo(datahist)
        
        return w.var('mean').getVal(),w.var('mean').getError(),w.var('sigma').getVal(),w.var('sigma').getError()
    
    
    


