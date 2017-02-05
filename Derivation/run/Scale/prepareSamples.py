import ROOT
ROOT.gSystem.Load("libKaMuCaDerivation")

Era='13TEV_Rereco'


#processor_UPSILON_data = ROOT.NtupleProcessor("UpsilonInput.root",1,9.452,0.0,0.0);
#processor_UPSILON_data.processTree("/scratch3/Kalman/data/{era}/Upsilon.root".format(era=Era),"pt1>4.0&&pt2>4.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>9.2&&mass<9.6")
#processor_UPSILON_data.close()

processor_JPSI_data = ROOT.NtupleProcessor("JPsiInput.root",1,3.09268,0.0,0.0,"DATA_Moriond17_13TeV");
processor_JPSI_data.processTree("/scratch3/Kalman/data/{era}/Jpsi.root".format(era=Era),"pt1>6&&pt2>6&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>2.895&&mass<3.295")
processor_JPSI_data.close()

processor_Z_data = ROOT.NtupleProcessor("ZInput.root",1,90.952,-7.929e-4,2.5,"DATA_Moriond17_13TeV");
processor_Z_data.processTree("/scratch3/Kalman/data/{era}/Z.root".format(era=Era),"pt1>20.0&&pt2>20.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>86&&mass<96")
processor_Z_data.close()



processor_JPSI_MC = ROOT.NtupleProcessor("JPsiMCInput.root",0,3.09268,0.0,0.0,"MC_Moriond17_13TeV");
processor_JPSI_MC.processTree("/scratch3/Kalman/data/{era}/JPsiMC.root".format(era=Era),"pt1>6&&pt2>6&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0&&mass>2.895&&mass<3.295")
processor_JPSI_MC.close()

processor_JPSI_MC = ROOT.NtupleProcessor("JPsiMCInput2.root",0,3.09268,0.0,0.0,"MC_Moriond17_13TeV",1);
processor_JPSI_MC.processTree("/scratch3/Kalman/data/{era}/JPsiMC20.root".format(era=Era),"pt1>6.0&&pt2>6.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0&&mass>2.895&&mass<3.295")
processor_JPSI_MC.close()

processor_Z_MC = ROOT.NtupleProcessor("ZMCInput.root",0,90.952,-7.929e-4,2.5,"MC_Moriond17_13TeV");
processor_Z_MC.processTree("/scratch3/Kalman/data/{era}/ZMC.root".format(era=Era),"pt1>20.0&&pt2>20.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0&&mass>86&&mass<96")
processor_Z_MC.close()

processor_Z_MC = ROOT.NtupleProcessor("ZMCInput2.root",0,90.952,-7.929e-4,2.5,"MC_Moriond17_13TeV");
processor_Z_MC.processTree("/scratch3/Kalman/data/{era}/ZJETSMC.root".format(era=Era),"pt1>20.0&&pt2>20.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0&&mass>86&&mass<96")
processor_Z_MC.close()












