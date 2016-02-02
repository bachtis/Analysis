import ROOT
import math
from KaMuCa.Derivation.tools.Partitioner import Partitioner
from KaMuCa.Derivation.tools.DataSetSplitter import DataSetSplitter

#Geometry map
pmap = Partitioner([1.0/150.,1.0/3.0],[-2.5,2.5],[-math.pi,math.pi],"")

def fx(line):
    return 1.0/10.0

def fy(line):
    return 0.0

def fz(line):
    return 0.0


#selectors
def jSelect(line):
    m = line.find("mass").getVal()
    mErr = line.find("massErr").getVal()

    pt1=1.0/line.find("c1").getVal()
    e1=line.find("eta1").getVal()
    p1=line.find("phi1").getVal()
    err1=line.find("cErr1").getVal()
    
    pt2=1.0/line.find("c2").getVal()
    e2=line.find("eta2").getVal()
    p2=line.find("phi2").getVal()
    err2=line.find("cErr2").getVal()


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


    if m>3.0 and m<3.2: 
        return True
    else:
        return False

def ySelect(line):
    m = line.find("mass").getVal()
    mErr = line.find("massErr").getVal()


    if m>9.2 and m<9.7 and mErr>0:
        return True
    else:
        return False



#Load Modifier that calibrates for the magnetic field
from KaMuCa.Derivation.tools.Calibrator import Calibrator
calibratorData = Calibrator("DATA_Prompt_13TeV",False,False,True)
splitter = DataSetSplitter(pmap,fx,fy,fz)
#splitter.split('/data/bachtis/CALIB/13TeV_76X/JDATA.root','JDATA13TeVHigh.root',[calibratorData.correct],[jSelect])
splitter.split('/data/bachtis/CALIB/13TeV_76X/JMCLOW.root','JMCLOW.root',[],[jSelect])







