import ROOT
ROOT.gSystem.Load("libKaMuCaDerivation")

Era='13TeV_2017'


#processor_UPSILON_data = ROOT.NtupleProcessor("UpsilonInput.root",1,9.452,0.0,0.0);
#processor_UPSILON_data.processTree("/scratch3/Kalman/data/{era}/Upsilon.root".format(era=Era),"pt1>4.0&&pt2>4.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>9.2&&mass<9.6")
#processor_UPSILON_data.close()


processor_Z_data = ROOT.NtupleProcessor("ZInput.root",1,90.851,0.0,2.5,"DATA_Moriond17_13TeV");
processor_Z_data.processTree("/scratch3/Kalman/data/{era}/Z.root".format(era=Era),"pt1>22.0&&pt2>22.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&mass>86&&mass<96&&cErr1<0.035&&cErr1<0.035")
processor_Z_data.close()

#processor_JPSI_MC = ROOT.NtupleProcessor("JPsiMCInputPre.root",0,3.095*0.9992,0.0,0.0,"MC_Moriond17_13TeV",3);
#processor_JPSI_MC.processTree("/scratch3/Kalman/data/{era}/JPsiMC.root".format(era=Era),"pt1>5&&pt2>5&&gc1>0&&gc2>0&&mass>2.895&&mass<3.295",2.4)
#processor_JPSI_MC.close()


#processor_JPSI_MC = ROOT.NtupleProcessor("JPsiInputRaw.root",1,3.095*0.9992,0.0,0.0,"DATA_Moriond17_13TeV",0);
#processor_JPSI_MC.processTree("/scratch3/Kalman/data/{era}/Jpsi.root".format(era=Era),"pt1>5&&pt2>5&&mass>2.895&&mass<3.295",2.4)
#processor_JPSI_MC.close()

processor_JPSI = ROOT.NtupleProcessor("JPsiInput.root",1,3.095,0.0,0.0,"DATA_Moriond17_13TeV");
processor_JPSI.processTree("/scratch3/Kalman/data/{era}/JPsi.root".format(era=Era),"pt1>6&&pt2>6&&mass>2.895&&mass<3.295",2.4)
processor_JPSI.close()





#processor_Z_MC = ROOT.NtupleProcessor("ZMCInput.root",0,90.851,0,2.5,"MC_Moriond17_13TeV");
#processor_Z_MC.processTree("/scratch3/Kalman/data/{era}/ZMC.root".format(era=Era),"pt1>22.0&&pt2>22.0&&abs(eta1)<1.6&&abs(eta2)<1.6&&gc1>0&&gc2>0&&mass>81&&mass<101&&cErr1<0.035&&cErr2<0.035",2.4)
#processor_Z_MC.close()




#1 For precalibration
#processor_Z_MC = ROOT.NtupleProcessor("ZMCInput.root",0,91.186/1.0027,0,2.5,"MC_Moriond17_13TeV",0);
#processor_Z_MC.processTree("/scratch3/Kalman/data/{era}/ZMC.root".format(era=Era),"pt1>22.0&&pt2>22.0&&abs(eta1)<2.4&&abs(eta2)<2.4&&gc1>0&&gc2>0",2.5)
#processor_Z_MC.close()










