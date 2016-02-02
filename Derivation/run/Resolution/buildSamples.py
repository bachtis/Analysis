import ROOT
import math
from defs import *
from KaMuCa.Derivation.tools.DataSetSplitter import DataSetSplitter


def fx(line):
    pt1 = 1.0/line.find('c1').getVal()
    return pt1*pt1

def fy(line):
    return (line.find('eta1').getVal())

def fz(line):
    return 0.0





#from KaMuCa.Derivation.tools.Calibrator import Calibrator
#calibratorData = Calibrator("DATA_Prompt_13TeV",False,False,True)


class etaselector(object):
    def __init__(self,pmap):
        self.pmap=pmap


    def select(self,line):
        pt1=1.0/line.find("c1").getVal()
        e1=line.find("eta1").getVal()
        p1=line.find("phi1").getVal()
        err1=line.find("cErr1").getVal()
       
        pt2=1.0/line.find("c2").getVal()
        e2=line.find("eta2").getVal()
        p2=line.find("phi2").getVal()
        err2=line.find("cErr2").getVal()





        bin1=self.pmap.binFromVals(pt1*pt1,e1,p1)
        bin2=self.pmap.binFromVals(pt2*pt2,e2,p2)

        x,be1,y = pmap.binXYZ(bin1)
        x,be2,y = pmap.binXYZ(bin2)

        err=0.5*math.sqrt(err1*err1+err2*err2)

        if err>0.2 or err==0.0:
            return False
        if be1==be2:
            return True
        return False

sel=etaselector(pmap)
from KaMuCa.Derivation.tools.Calibrator import Calibrator
calibratorData = Calibrator("DATA_76X_13TeV",True,False,False)


splitter = DataSetSplitter(pmap,fx,fy,fz)   
splitter.split('/data/bachtis/CALIB/13TeV_76X/JDATA.root','JDATA.root',[calibratorData.correct],[sel.select])
splitter = DataSetSplitter(pmap,fx,fy,fz)   
splitter.split('/data/bachtis/CALIB/13TeV_76X/ZDATA.root','ZDATA.root',[calibratorData.correct],[sel.select])









