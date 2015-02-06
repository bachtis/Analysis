import ROOT

def fetch(filename,histo,name,f):
    ff=ROOT.TFile(filename)
    hh=ff.Get(histo).Clone()
    hh.SetName(name)
    f.cd()
    hh.Write()


f=ROOT.TFile("dataInputs.root","RECREATE")
fetch("step0_magneticMap/mapCalibration.root","mapCorrection","magnetic",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","A_0","A1",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","K_0","A2",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","B_0","B",f)
fetch("step3_KalmanFilter/kalmanScale_data.root","M_0","e",f)
fetch("Resolution/fitResults.root","A_data","sigma_A_target",f)
fetch("Resolution/fitResults.root","B_data","sigma_B_target",f)
fetch("Resolution/fitResults.root","C_data","sigma_C_target",f)
fetch("Resolution/fitResults.root","A_mc","sigma_A_ref",f)
fetch("Resolution/fitResults.root","B_mc","sigma_B_ref",f)
fetch("Resolution/fitResults.root","C_mc","sigma_C_ref",f)
fetch("Resolution/fitResultsEbE.root","A_data","ebe_A",f)
fetch("Resolution/fitResultsEbE.root","B_data","ebe_B",f)
fetch("Resolution/fitResultsEbE.root","C_data","ebe_C",f)

f.Close()

f=ROOT.TFile("mcInputs.root","RECREATE")
fetch("step3_KalmanFilter/kalmanScale_mc.root","A_0","A1",f)
fetch("step3_KalmanFilter/kalmanScale_mc.root","K_0","A2",f)
fetch("step3_KalmanFilter/kalmanScale_mc.root","B_0","B",f)
fetch("step3_KalmanFilter/kalmanScale_mc.root","M_0","e",f)
fetch("Resolution/fitResults.root","A_data","sigma_A_ref",f)
fetch("Resolution/fitResults.root","B_data","sigma_B_ref",f)
fetch("Resolution/fitResults.root","C_data","sigma_C_ref",f)
fetch("Resolution/fitResults.root","A_mc","sigma_A_target",f)
fetch("Resolution/fitResults.root","B_mc","sigma_B_target",f)
fetch("Resolution/fitResults.root","C_mc","sigma_C_target",f)
fetch("Resolution/fitResultsEbE.root","A_mc","ebe_A",f)
fetch("Resolution/fitResultsEbE.root","B_mc","ebe_B",f)
fetch("Resolution/fitResultsEbE.root","C_mc","ebe_C",f)

f.Close()










