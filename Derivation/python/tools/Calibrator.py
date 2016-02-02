import ROOT
import math

class Calibrator(object):
    def __init__(self,corr,ptCalib=False,errCalib=False,bOnlyCalib=False,msOnlyCalib=False):
        self.corr=corr
        self.ptCalib = ptCalib
        self.errCalib = errCalib
        self.bOnlyCalib = bOnlyCalib
        self.calibrator = ROOT.KalmanMuonCalibrator(corr)
        self.msOnlyCalib=msOnlyCalib

    def correct(self,line):
        pt1 = 1./line.find('c1').getVal()
        eta1 = line.find('eta1').getVal()
        phi1 = line.find('phi1').getVal()
        err1=line.find('cErr1').getVal()
        pt2 = 1./line.find('c2').getVal()
        eta2 = line.find('eta2').getVal()
        phi2 = line.find('phi2').getVal()
        err2=line.find('cErr2').getVal()
        m=line.find('mass').getVal()

        
        if self.msOnlyCalib:
            err1 = self.calibrator.getCorrectedErrorMS(pt1,eta1,err1)
            err2 = self.calibrator.getCorrectedErrorMS(pt2,eta2,err2)
        else:    
            if self.errCalib:
                err1 = self.calibrator.getCorrectedError(pt1,eta1,err1)
                err2 = self.calibrator.getCorrectedError(pt2,eta2,err2)
             
        if self.bOnlyCalib:
            pt1 = self.calibrator.getCorrectedPtMag(pt1,eta1,phi1)
            pt2 = self.calibrator.getCorrectedPtMag(pt2,eta2,phi2)
        else:
            if self.ptCalib:
                pt1 = self.calibrator.getCorrectedPt(pt1,eta1,phi1,1)
                pt2 = self.calibrator.getCorrectedPt(pt2,eta2,phi2,-1)



            
        line.find('c1').setVal(1.0/pt1)
        line.find('c2').setVal(1.0/pt2)
        line.find('cErr1').setVal(err1)
        line.find('cErr2').setVal(err2)
        line.find('massErr').setVal(0.5*math.sqrt(err1*err1+err2*err2)*m)

        #recalculate the mass
        v1=ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(1./line.find('c1').getVal(),
                        line.find('eta1').getVal(),
                        line.find('phi1').getVal(),
                        0.1056583715)
        
        v2=ROOT.TLorentzVector()
        v2.SetPtEtaPhiM(1./line.find('c2').getVal(),
                        line.find('eta2').getVal(),
                        line.find('phi2').getVal(),
                        0.1056583715)
        line.find('mass').setVal((v1+v2).M())


    def smear(self,line):
        pt1 = 1./line.find('c1').getVal()
        eta1 = line.find('eta1').getVal()
        phi1 = line.find('phi1').getVal()
        err1=line.find('cErr1').getVal()
        pt2 = 1./line.find('c2').getVal()
        eta2 = line.find('eta2').getVal()
        phi2 = line.find('phi2').getVal()
        err2=line.find('cErr2').getVal()
        m=line.find('mass').getVal()


        
#        print 'before',pt1,pt2,eta1,eta2
        pt1 = self.calibrator.smear(pt1,eta1)
        pt2 = self.calibrator.smear(pt2,eta2)
#        print 'after',pt1,pt2
        
            
        line.find('c1').setVal(1.0/pt1)
        line.find('c2').setVal(1.0/pt2)


        #recalculate the mass
        v1=ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(1./line.find('c1').getVal(),
                        line.find('eta1').getVal(),
                        line.find('phi1').getVal(),
                        0.1056583715)
        
        v2=ROOT.TLorentzVector()
        v2.SetPtEtaPhiM(1./line.find('c2').getVal(),
                        line.find('eta2').getVal(),
                        line.find('phi2').getVal(),
                        0.1056583715)
        line.find('mass').setVal((v1+v2).M())


    def smearUsingEbE(self,line):
        pt1 = 1./line.find('c1').getVal()
        eta1 = line.find('eta1').getVal()
        phi1 = line.find('phi1').getVal()
        err1=line.find('cErr1').getVal()
        pt2 = 1./line.find('c2').getVal()
        eta2 = line.find('eta2').getVal()
        phi2 = line.find('phi2').getVal()
        err2=line.find('cErr2').getVal()
        m=line.find('mass').getVal()


        
#        print 'before',pt1,pt2,eta1,eta2
        pt1 = self.calibrator.smearUsingEbE(pt1,eta1,err1)
        pt2 = self.calibrator.smearUsingEbE(pt2,eta2,err2)
#        print 'after',pt1,pt2
        
            
        line.find('c1').setVal(1.0/pt1)
        line.find('c2').setVal(1.0/pt2)


        #recalculate the mass
        v1=ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(1./line.find('c1').getVal(),
                        line.find('eta1').getVal(),
                        line.find('phi1').getVal(),
                        0.1056583715)
        
        v2=ROOT.TLorentzVector()
        v2.SetPtEtaPhiM(1./line.find('c2').getVal(),
                        line.find('eta2').getVal(),
                        line.find('phi2').getVal(),
                        0.1056583715)
        line.find('mass').setVal((v1+v2).M())
        



