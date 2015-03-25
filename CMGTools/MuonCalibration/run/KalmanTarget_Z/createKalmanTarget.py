import ROOT
from CMGTools.MuonCalibration.tools.Smearing  import smearAbsolute,smearEbE2D

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()

ROOT.gSystem.Load("libCMGToolsMuonCalibration")


f=ROOT.TFile("../../data/ZGEN.root")
data=f.Get("data").reduce("massRaw>85.&&massRaw<95.&&abs(etaRaw1)<2.4&&abs(etaRaw2)<2.4")
#data=smearAbsolute(data,False)
data=smearEbE2D(data,1,1,True)
profile= ROOT.TProfile("mass_fit","MassFit",80,-2,2,85,95,"s")

for i in range(0,data.numEntries()):
    line=data.get(i)
    rapidity = line.find("rapidity").getVal()
    mass = line.find("massRaw").getVal()
    profile.Fill(rapidity,mass)


f2 = ROOT.TFile("kalmanTargetZ_Ext.root","RECREATE")
f2.cd()
profile.Write()
f2.Close()



