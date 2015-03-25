import ROOT
import math
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()

ROOT.gSystem.Load("libCMGToolsMuonCalibration")


from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.Unfolding import Unfolding,AnalyticalUnfolding1D,AnalyticalUnfolding,AnalyticalUnfoldingMass
from CMGTools.MuonCalibration.tools.calibTools import correctDataSet
from CMGTools.MuonCalibration.tools.Scheduler import Scheduler
from CMGTools.MuonCalibration.tools.LineshapeFitter import LineshapeFitter
from CMGTools.MuonCalibration.tools.Smearing import smearAbsolute,smearEbE2D,smearRelative


#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)


curvArrZ=[20.,30.,40.0,50.0,60.0,80.0]
curvArrJ=[3.0,5.0,7.0,10.0,12.5,15.0]
curvArrY=[3.0,5.0,7.0,10.0,12.5,15.0,20.0]
etaArr=[0.0,0.9,2.5]
phiArr=[-math.pi,math.pi]






