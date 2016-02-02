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
    pt1 = 1.0/line.find("c1").getVal()
    pt2 = 1.0/line.find("c2").getVal()


    if m>3.0 and m<3.2 and mErr>0 and pt1>3 and pt2>3 and pt1<50 and pt2<50:
        return True
    else:
        return False

#Load Modifier that calibrates for the magnetic field
from KaMuCa.Derivation.tools.Calibrator import Calibrator
calibratorData = Calibrator("DATA_53X_8TeV",False,True,True)
calibratorMC = Calibrator("MC_53X_8TeV",False,True,False)
splitter = DataSetSplitter(pmap,fx,fy,fz)

splitter.split('/data/bachtis/CALIB/8TeVNew/JDATA.root','JDATA8TeV.root',[calibratorData.correct],[jSelect])
splitter.split('/data/bachtis/CALIB/8TeVNew/JMC.root','JMC8TeV.root',[calibratorMC.correct],[jSelect])



#splitter.split('../../data/JMC.root','JMC.root',[],[jSelect])
#splitter.split('../../data/YMC.root','YMC.root',[],[ySelect])






