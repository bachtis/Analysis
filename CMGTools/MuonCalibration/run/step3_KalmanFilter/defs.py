import ROOT
import math
ROOT.gSystem.Load("libCMGToolsMuonCalibration")


from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.KalmanCalibrator  import *


#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

#w.var('massRaw').setBins(100)
#w.var('massRaw').setMin(2.8)
#w.var('massRaw').setMax(3.4)



curvArr= [1/100.,1/4.]

phiArr=[]
for i in range(0,11):
    phiArr.append(-math.pi+2*math.pi*i/10.)

etaArr=[]
for i in range(0,11):
    etaArr.append(-0.9+1.8*i/10.)



pmap = PartitionMap(curvArr,etaArr,phiArr,"")

#pmap.declareData('A',1.00)
#pmap.declareData('B',0.00)
#pmap.declareData('M',0.00)





