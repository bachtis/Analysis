import ROOT
import math
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libCMGToolsMuonCalibration")


from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.KalmanCalibrator  import *
from CMGTools.MuonCalibration.tools.Scheduler import Scheduler
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter
from CMGTools.MuonCalibration.tools.Smearing import smearAbsolute,smearEbE2D
from CMGTools.MuonCalibration.tools.calibTools import correctDataSet


#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

w.var('massRaw').setBins(100)
w.var('massRaw').setMin(2.8)
w.var('massRaw').setMax(3.4)



curvArr= [1./100, 1./5.5]
phiArr=[-math.pi,math.pi]

etaArr=[]
for i in range(0,21):
    etaArr.append(-0.9+1.8*i/20)


pmap = PartitionMap(curvArr,etaArr,phiArr,"")





