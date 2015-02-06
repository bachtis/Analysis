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
from CMGTools.MuonCalibration.tools.Smearing import smearAbsolute,smearEbE2D,smearFlat

#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

w.var('massRaw').setBins(100)
w.var('massRaw').setMin(80)
w.var('massRaw').setMax(120)


curvArr=[1./100. ,1./50. , 1./35. ,1./25, 1./16. , 1./14. , 1./12. , 1./10. , 1./8., 1./6.5, 1/5.]
phiArr=[-math.pi,math.pi]
etaArr = [-1.1,-0.9,-0.5,0.0,0.5,0.9,1.1]
pmap = PartitionMap(curvArr,etaArr,phiArr,"")



