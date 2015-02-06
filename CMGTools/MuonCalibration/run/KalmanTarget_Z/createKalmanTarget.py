import ROOT
from CMGTools.MuonCalibration.tools.Smearing  import smearAbsolute

ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.AutoLibraryLoader.enable()

ROOT.gSystem.Load("libCMGToolsMuonCalibration")


f=ROOT.TFile("../../data/ZGEN.root")
data=f.Get("data").reduce("massRaw>85&&massRaw<95&&abs(etaRaw1)<0.9&&abs(etaRaw2)<0.9")
data=smearAbsolute(data,False)
profile= ROOT.TProfile("mass_fit","MassFit",50,-0.9,0.9,85,95,"s")


for i in range(0,data.numEntries()):
    line=data.get(i)
    rapidity = line.find("rapidity").getVal()
    mass = line.find("massRaw").getVal()
    profile.Fill(rapidity,mass)


f2 = ROOT.TFile("kalmanTargetZ_MC.root","RECREATE")
f2.cd()
profile.Write()
f2.Close()



