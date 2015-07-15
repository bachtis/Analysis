import ROOT
import math

def smearAbsolute(data,isData = True):
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

        pt1 = calibrator.smearGEN(pt1,eta1)
        pt2 = calibrator.smearGEN(pt2,eta2)

            
        line.find('curvRaw1').setVal(1.0/pt1)
        line.find('curvRaw2').setVal(1.0/pt2)

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


def smearRelative(data,isData):
    calibrator = ROOT.KalmanCalibratorParam(isData)
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
        print 'before',pt1,pt2
        pt1 = calibrator.smear(pt1,eta1)
        pt2 = calibrator.smear(pt2,eta2)
        print 'after',pt1,pt2

            
        line.find('curvRaw1').setVal(1.0/pt1)
        line.find('curvRaw2').setVal(1.0/pt2)

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



def smearEbE2D(data,shift = 1.0,errorScale=1.0,updateMass = True):
    random = ROOT.TRandom3(101082)
    newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
    for evt in range(0,data.numEntries()):
        line = data.get(evt)
       
        errx = errorScale*2*line.find('massErrRaw1').getVal()*line.find('curvRaw1').getVal()/line.find('massRaw').getVal()
        erry = errorScale*2*line.find('massErrRaw2').getVal()*line.find('curvRaw2').getVal()/line.find('massRaw').getVal()
        
        s1=random.Gaus(0.0,errx)
        s2=random.Gaus(0.0,erry)

        c1 = shift*line.find('curvRaw1').getVal()+s1
        c2 = shift*line.find('curvRaw2').getVal()+s2

        if c1<=0:
            continue
        if c2<=0:
            continue
        
        line.find('curvRaw1').setVal(c1)
        line.find('curvRaw2').setVal(c2)

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

        if updateMass:
            line.find('massRaw').setVal((v1+v2).M())

            


        newData.add(line)

    return newData


def smearFlat(data,shift = 1.0,error=1.0,updateMass = True):
    random = ROOT.TRandom3(101082)
    newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
    for evt in range(0,data.numEntries()):
        line = data.get(evt)
       
        s1=random.Gaus(0.0,error)
        s2=random.Gaus(0.0,error)

        c1 = shift*line.find('curvRaw1').getVal()+s1
        c2 = shift*line.find('curvRaw2').getVal()+s2

        if c1<=0:
            continue
        if c2<=0:
            continue
        
        line.find('curvRaw1').setVal(c1)
        line.find('curvRaw2').setVal(c2)

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

        if updateMass:
            line.find('massRaw').setVal((v1+v2).M())

            


        newData.add(line)

    return newData




def createDataLike(data,gen):
    random = ROOT.TRandom3(101082)
    newData = ROOT.RooDataSet(data.GetName()+'cal',data.GetName(),data.get())
    genEntries = gen.numEntries() 

    for evt in range(0,data.numEntries()):
        line = data.get(evt)

        genLine = gen.get(int(random.Rndm()*genEntries))
        line.find('massRaw').setVal(genLine.find('massRaw').getVal())
        newData.add(line)

    return newData

