import ROOT
import math
def correctDataSet(data,isData=True,ptCalib=False,errCalib=False,bOnlyCalib=False):
    calibrator = ROOT.KalmanCalibrator(isData)
    newData = ROOT.RooDataSet("data","data",data.get())
    for i in range(0,data.numEntries()):
        line =data.get(i)
        pt1 = 1./line.find('curvRaw1').getVal()
        eta1 = line.find('etaRaw1').getVal()
        phi1 = line.find('phiRaw1').getVal()
        err1=line.find('massErrRaw1').getVal()
        pt2 = 1./line.find('curvRaw2').getVal()
        eta2 = line.find('etaRaw2').getVal()
        phi2 = line.find('phiRaw2').getVal()
        err2=line.find('massErrRaw2').getVal()
        m=line.find('massRaw').getVal()


        if errCalib:    
            err1 = calibrator.getCorrectedError(pt1,eta1,err1/m)*m
            err2 = calibrator.getCorrectedError(pt2,eta2,err2/m)*m
       

             
        if bOnlyCalib:
            pt1 = calibrator.getCorrectedPtMag(pt1,eta1,phi1)
            pt2 = calibrator.getCorrectedPtMag(pt2,eta2,phi2)
        else:
            if ptCalib:
                pt1 = calibrator.getCorrectedPt(pt1,eta1,phi1,1)
                pt2 = calibrator.getCorrectedPt(pt2,eta2,phi2,-1)



            
        line.find('curvRaw1').setVal(1.0/pt1)
        line.find('curvRaw2').setVal(1.0/pt2)
        line.find('massErrRaw1').setVal(err1)
        line.find('massErrRaw2').setVal(err2)
        line.find('massErrRaw').setVal(math.sqrt(err1*err1+err2*err2))

        #recalculate the mass
        v1=ROOT.TLorentzVector()
        v1.SetPtEtaPhiM(1./line.find('curvRaw1').getVal(),
                        line.find('etaRaw1').getVal(),
                        line.find('phiRaw1').getVal(),
                        0.1056583715)
        
        v2=ROOT.TLorentzVector()
        v2.SetPtEtaPhiM(1./line.find('curvRaw2').getVal(),
                        line.find('etaRaw2').getVal(),
                        line.find('phiRaw2').getVal(),
                        0.1056583715)
        line.find('massRaw').setVal((v1+v2).M())
        
        newData.add(line)
            
    return newData

