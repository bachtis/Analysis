import ROOT

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()
ROOT.gSystem.Load("libCMGToolsMuonCalibration")
import math
from CMGTools.MuonCalibration.tools.calibTools import correctDataSetParam
from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter
from CMGTools.MuonCalibration.tools.EbEFitter import EbEFitter
from CMGTools.MuonCalibration.tools.Scheduler import Scheduler
from CMGTools.MuonCalibration.tools.Smearing import smearAbsolute,smearEbE2D,smearFlat

#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

w.var('massRaw').setBins(100)
w.var('massRaw').setMin(80)
w.var('massRaw').setMax(120)


curvArr=[3.0,5.0,7.0,10.0,20.0,30.0,40.0,50.0,100.0]
etaArr=[0.0,0.25,0.5,0.7,0.9,1.1,1.4,1.7,2.0,2.4]
phiArr=[-math.pi,math.pi]
pmap = PartitionMap(curvArr,etaArr,phiArr,"")


