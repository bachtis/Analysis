import ROOT
ROOT.gSystem.Load("libKaMuCaDerivation")

Era='13TEV_Rereco'

processor_JPSI_data = ROOT.NtupleProcessor("JPsiInput.root",1,3.095,0.0,0.0);
processor_JPSI_data.processTree("/scratch3/Kalman/data/{era}/Jpsi.root".format(era=Era),"pt1>3.0&&pt2>3.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>2.895&&mass<3.295")
processor_JPSI_data.Write()


processor_UPSILON_data = ROOT.NtupleProcessor("UpsilonInput.root",1,9.452,0.0,0.0);
processor_UPSILON_data.processTree("/scratch3/Kalman/data/{era}/Upsilon.root".format(era=Era),"pt1>3.0&&pt2>3.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>9.2&&mass<9.6")
processor_UPSILON_data.Write()






