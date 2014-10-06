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
w.var('massRaw').setMin(70)
w.var('massRaw').setMax(120)



curvArr= [0.0,9999999.0]

phiArr=[-math.pi,math.pi]

etaArr=[]
for i in range(0,51):
    etaArr.append(-5+10*i/50.)

pmap = PartitionMap(curvArr,etaArr,phiArr,"results.root")







