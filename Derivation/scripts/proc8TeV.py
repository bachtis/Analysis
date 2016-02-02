

import ROOT 




p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/8TeVNew/JMC.root",False,False)
p.processTree("/tmp/bachtis/JMC.root","MuPos_pt>3&&MuNeg_pt>3&&(MuPos_charge+MuNeg_charge)==0&&Z_mass>2.9&&Z_mass<3.3&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02")
p.write()



p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/8TeVNew/JGEN.root",False,True)
p.processTree("/tmp/bachtis/JMC.root","MuPos_pt>3&&MuNeg_pt>3&&(MuPos_charge+MuNeg_charge)==0&&Z_mass>2.9&&Z_mass<3.3&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02&&MuPosGenStatus1_pt>0&&MuNegGenStatus1_pt>0")
p.write()


p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/8TeVNew/JDATA.root",True,False)
p.processTree("/tmp/bachtis/JDATA.root","MuPos_pt>3&&MuNeg_pt>3&&(MuPos_charge+MuNeg_charge)==0&&Z_mass>2.9&&Z_mass<3.3&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02")
p.write()



p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/8TeVNew/ZMC.root",False,False)
p.processTree("/tmp/bachtis/ZMC.root","MuPos_pt>20&&MuNeg_pt>20&&(MuPos_charge+MuNeg_charge)==0&&Z_mass>80&&Z_mass<110&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/8TeVNew/ZGEN.root",False,True)
p.processTree("/tmp/bachtis/ZMC.root","MuPos_pt>20&&MuNeg_pt>20&&(MuPos_charge+MuNeg_charge)==0&&Z_mass>80&&Z_mass<110&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02&&MuPosGenStatus1_pt>0&&MuNegGenStatus1_pt>0")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/8TeVNew/ZDATA.root",True,False)
p.processTree("/tmp/bachtis/ZDATA.root","MuPos_pt>20&&MuNeg_pt>20&&(MuPos_charge+MuNeg_charge)==0&&Z_mass>80&&Z_mass<110&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02")
p.write()


