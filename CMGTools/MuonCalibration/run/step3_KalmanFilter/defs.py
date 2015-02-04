import ROOT
import math

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()

ROOT.gSystem.Load("libCMGToolsMuonCalibration")

from CMGTools.MuonCalibration.tools.calibTools import correctDataSet
from CMGTools.MuonCalibration.tools.workspaceTools import prepareWorkspace
from CMGTools.MuonCalibration.tools.PartitionMap import PartitionMap
from CMGTools.MuonCalibration.tools.DataSetBuilder import DataSetBuilder
from CMGTools.MuonCalibration.tools.Unfolding import AnalyticalUnfolding,Unfolding
from CMGTools.MuonCalibration.tools.Smearing import smearAbsolute,smearEbE2D,smearFlat


#create a MC dataset for the lineshape
w=ROOT.RooWorkspace('w','w')
prepareWorkspace(w)

curvArr= [1/100.,1/5.5]

phiArr=[]
for i in range(0,11):
    phiArr.append(-math.pi+2*math.pi*i/10.)

etaArr=[]
for i in range(0,11):
    etaArr.append(-0.9+1.8*i/10.)




pmap = PartitionMap(curvArr,etaArr,phiArr,"")




