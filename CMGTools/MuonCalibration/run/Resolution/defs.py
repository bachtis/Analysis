import ROOT

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libCMGToolsMuonCalibration")
import math
from CMGTools.MuonCalibration.tools.calibTools import correctDataSet
from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter
from CMGTools.MuonCalibration.tools.Scheduler import Scheduler

#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

w.var('massRaw').setBins(100)
w.var('massRaw').setMin(80)
w.var('massRaw').setMax(120)


curvArr=[1/100.,1/5.]
phiArr=[-math.pi,math.pi]
etaArr = []
for i in range(0,21):
    etaArr.append(-1.1+2.2*i/20.)
pmap = PartitionMap(curvArr,etaArr,phiArr,"")







