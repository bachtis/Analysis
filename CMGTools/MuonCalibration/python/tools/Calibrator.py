import ROOT


class Calibrator:
    def __init__(self,filename,iteration):
        self.f = ROOT.TFile(filename)
        print 'Cpos_{i}'.format(i=iteration)
        self.pos = self.f.Get('Cpos_{i}'.format(i=iteration))
        self.neg = self.f.Get('Cneg_{i}'.format(i=iteration))

    def calibrateP(self,charge,pt,eta,phi):
        curv = 1/pt
        if charge >0:
            binx=self.pos.GetXaxis().FindBin(curv)
            biny=self.pos.GetYaxis().FindBin(eta)
            binz=self.pos.GetZaxis().FindBin(phi)
            bin = self.pos.GetBin(binx,biny,binz)

            factor = self.pos.GetBinContent(bin)
            if factor>0:
                return pt/factor, eta,phi,curv*factor
            else:
                return pt, eta,phi,curv

        else:
            binx=self.neg.GetXaxis().FindBin(curv)
            biny=self.neg.GetYaxis().FindBin(eta)
            binz=self.neg.GetZaxis().FindBin(phi)
            bin = self.neg.GetBin(binx,biny,binz)
            factor = self.neg.GetBinContent(bin)
            if factor>0:
                return pt/factor, eta,phi,curv*factor
            else:
                return pt, eta,phi,curv


    def calibrateDataset(self,data,w):
        newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),ROOT.RooArgSet(
            w.var('curvRaw1'),
            w.var('etaRaw1'),
            w.var('phiRaw1'),
            w.var('curvRaw2'),
            w.var('etaRaw2'),
            w.var('phiRaw2'),
            w.var('massRaw')
            ))
        
        for evt in range(0,data.numEntries()):
            argset = data.get(evt)

            pt1 = 1./argset.find('curvRaw1').getVal()
            eta1 = argset.find('etaRaw1').getVal()
            phi1 = argset.find('phiRaw1').getVal()

            pt1,eta1,phi1,curv1 = self.calibrateP(1,pt1,eta1,phi1)
            
            
            pt2 = 1./argset.find('curvRaw2').getVal()
            eta2 = argset.find('etaRaw2').getVal()
            phi2 = argset.find('phiRaw2').getVal()


            pt2,eta2,phi2,curv2 = self.calibrateP(-1,pt2,eta2,phi2)
            

            w.var('curvRaw1').setVal(1./pt1)    
            w.var('curvRaw2').setVal(1./pt2)    
            w.var('etaRaw1').setVal(eta1)    
            w.var('etaRaw2').setVal(eta2)    
            w.var('phiRaw1').setVal(phi1)    
            w.var('phiRaw2').setVal(phi2)    
                

            v1=ROOT.TLorentzVector()
            v2=ROOT.TLorentzVector()
            v1.SetPtEtaPhiM(pt1,eta1,phi1,w.var('muMass').getVal())
            v2.SetPtEtaPhiM(pt2,eta2,phi2,w.var('muMass').getVal())
            w.var('massRaw').setVal((v1+v2).M())
            newData.add(ROOT.RooArgSet(w.var('curvRaw1'),
                                       w.var('etaRaw1'),
                                       w.var('phiRaw1'),
                                       w.var('curvRaw2'),
                                       w.var('etaRaw2'),
                                       w.var('phiRaw2'),
                                       w.var('massRaw')
                                       ))
        return newData
        
