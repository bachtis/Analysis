import ROOT

def fetch(filename,histo,name,f):
    ff=ROOT.TFile(filename)
    hh=ff.Get(histo).Clone()
    hh.SetName(name)
    f.cd()
    hh.Write()


f=ROOT.TFile("dataInputs.root","RECREATE")
fetch("step0_magneticMap/mapCalibration.root","mapCorrection","magnetic",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","A_0","A1",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","K_0","A2",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","B_0","B",f)
fetch("KalmanFilter_Ext/kalmanScale_data.root","M_0","e",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_A_target","sigma_A_target",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_B_target","sigma_B_target",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_C_target","sigma_C_target",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_A_target","sigma_A_ref",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_B_target","sigma_B_ref",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","sigma_C_target","sigma_C_ref",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","ebe_A","ebe_A",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","ebe_B","ebe_B",f)
fetch("../../KalmanCalibration/data/dataInputsBKP.root","ebe_C","ebe_C",f)


f.Close()

f=ROOT.TFile("mcInputs.root","RECREATE")
fetch("KalmanFilter_Ext/kalmanScale_mc.root","A_0","A1",f)
fetch("KalmanFilter_Ext/kalmanScale_mc.root","K_0","A2",f)
fetch("KalmanFilter_Ext/kalmanScale_mc.root","B_0","B",f)
fetch("KalmanFilter_Ext/kalmanScale_mc.root","M_0","e",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_A_target","sigma_A_target",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_B_target","sigma_B_target",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_C_target","sigma_C_target",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_A_target","sigma_A_ref",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_B_target","sigma_B_ref",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","sigma_C_target","sigma_C_ref",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","ebe_A","ebe_A",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","ebe_B","ebe_B",f)
fetch("../../KalmanCalibration/data/mcInputsBKP.root","ebe_C","ebe_C",f)


f.Close()










