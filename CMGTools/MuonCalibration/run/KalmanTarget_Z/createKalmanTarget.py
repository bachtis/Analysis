import ROOT


f=ROOT.TFile("../../data/ZGEN.root")
data=f.Get("data").reduce("massRaw>85&&massRaw<95&&abs(etaRaw1)<0.9&&abs(etaRaw2)<0.9")
profile= ROOT.TProfile("mass_fit","MassFit",25,-0.9,0.9,80,100,"s")


for i in range(0,data.numEntries()):
    line=data.get(i)
    rapidity = line.find("rapidity").getVal()
    mass = line.find("massRaw").getVal()
    profile.Fill(rapidity,mass)


f2 = ROOT.TFile("kalmanTargetZ.root","RECREATE")
f2.cd()
profile.Write()
f2.Close()



