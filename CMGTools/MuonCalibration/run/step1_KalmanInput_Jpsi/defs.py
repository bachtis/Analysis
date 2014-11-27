import ROOT
import math
ROOT.gSystem.Load("libCMGToolsMuonCalibration")


from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.KalmanCalibrator  import *
from CMGTools.MuonCalibration.tools.Scheduler import Scheduler
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter


#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

w.var('massRaw').setBins(100)
w.var('massRaw').setMin(2.8)
w.var('massRaw').setMax(3.4)



curvArr= [0.0,999999999.]

phiArr=[-math.pi,math.pi]

etaArr=[]
for i in range(0,26):
    etaArr.append(-1.3+2.6*i/25.0)



pmap = PartitionMap(curvArr,etaArr,phiArr,"")





