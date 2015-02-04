import ROOT

def fetch(filename,histo,name,f):
    ff=ROOT.TFile(filename)
    hh=ff.Get(histo).Clone()
    hh.SetName(name)
    f.cd()
    hh.Write()


f=ROOT.TFile("dataInputs.root","RECREATE")
fetch("step0_magneticMap/mapCalibration.root","mapCorrection","magnetic",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","A_1","A1",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","K_1","A2",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","B_1","B",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","M_1","e",f)
fetch("Resolution/resolutionData.root","b2_0","b2",f)
fetch("Resolution/resolutionHighPtData.root","a2_0","a2",f)
fetch("Resolution/resolutionMaps.root","dataMap","resolution",f)
fetch("Resolution/resolutionMaps.root","diff","smearing",f)
f.Close()

f=ROOT.TFile("mcInputs.root","RECREATE")
fetch("step3_KalmanFilter/kalmanScale_mc.root","A_0","A1",f)
fetch("step3_KalmanFilter/kalmanScale_mc.root","K_0","A2",f)
fetch("step3_KalmanFilter/kalmanScale_mc.root","B_0","B",f)
fetch("step3_KalmanFilter/kalmanScale_mc.root","M_0","e",f)
fetch("Resolution/resolutionMC.root","b2_0","b2",f)
fetch("Resolution/resolutionHighPtMC.root","a2_0","a2",f)
fetch("Resolution/resolutionMaps.root","dataMap","resolution",f)
fetch("Resolution/resolutionMaps.root","diff","smearing",f)

f.Close()










