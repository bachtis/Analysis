import ROOT
import math

random = ROOT.TRandom3(101082)

f=ROOT.TFile("/data/bachtis/CALIB/8TeVNew/JGEN.root")
data=f.Get("data").reduce("abs(eta1)<2.5&&abs(eta2)<2.5")

LO=3.0
HI=3.2

from KaMuCa.Derivation.tools.Calibrator import Calibrator


profile= ROOT.TProfile("mass_fit","MassFit",46,-2.5,2.5,LO,HI)
histo2D= ROOT.TH2D("histo","MassFit",46,-2.5,2.5,20,LO,HI)
profileE= ROOT.TProfile("mass_fitE","MassFit",50,0,3e-4,LO,HI)
histo2DE= ROOT.TH2D("histoE","MassFit",50,0,3e-4,20,LO,HI)

profile2DE= ROOT.TProfile2D("mass_fit2D","MassFit",20,0,9e-4,23,-2.3,2.3,LO,HI)
histo3DE= ROOT.TH3F("histo3D","MassFit",20,0,9e-4,23,-2.3,2.3,20,LO,HI)

calibrator = Calibrator("MC_53X_8TeV",False,True,False)




random=ROOT.TRandom3(101082)


newData = ROOT.RooDataSet("data","data",data.get())
for i in range(0,data.numEntries()):
    line=data.get(i)
    calibrator.correct(line)

    v1=ROOT.TLorentzVector()
    v1.SetPtEtaPhiM(1./line.find('c1').getVal(),
                    line.find('eta1').getVal(),
                    line.find('phi1').getVal(),
                    0.1056583715)
    
    v2=ROOT.TLorentzVector()
    v2.SetPtEtaPhiM(1./line.find('c2').getVal(),
                    line.find('eta2').getVal(),
                    line.find('phi2').getVal(),
                    0.1056583715)
    
    v = v1+v2


    rapidity = v.Rapidity()
    massR = v.M()


    for j in range(1,profile2DE.GetXaxis().GetNbins()+1):
        err2 = profile2DE.GetXaxis().GetBinCenter(j)
        err = math.sqrt(err2)*massR
        mass=massR+random.Gaus(0.0,err)
        if mass>LO and mass<HI:
            profile2DE.Fill(err2,rapidity,mass)
            histo3DE.Fill(err2,rapidity,mass)


    err = line.find("massErr").getVal()
    if err<0:
        continue
    mass=massR+random.Gaus(0.0,err)
    if mass>LO and mass<HI:
        profile.Fill(rapidity,mass)
        histo2D.Fill(rapidity,mass)
        profileE.Fill(err*err/(mass*mass),mass)
        histo2DE.Fill(err*err/(mass*mass),mass)
        
        newData.add(line)


f2 = ROOT.TFile("kalmanTargetJ.root","RECREATE")
f2.cd()
profile.Write()
histo2D.Write()
profile2DE.Write()

profileE.Write()
histo2DE.Write()

histo3DE.Write()

newData.Write()
f2.Close()



