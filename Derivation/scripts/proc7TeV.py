

import ROOT 




p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/JMC.root",0,0)
p.processTree("/tmp/bachtis/JMC.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>2.9&&Z_mass<3.3&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/JGEN.root",0,1)
p.processTree("/tmp/bachtis/JMC.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>2.9&&Z_mass<3.3&&MuPosGenStatus1_pt>0&&MuNegGenStatus1_pt>0&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()


p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/JDATA.root",1,0)
p.processTree("/tmp/bachtis/JDATA.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>2.9&&Z_mass<3.3&&abs(MuPos_dxy)<0.02&&abs(MuNeg_dxy)<0.02&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()



p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/ZMC.root",0,0)
p.processTree("/tmp/bachtis/ZMC.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>80&&Z_mass<120&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/ZGEN.root",0,1)
p.processTree("/tmp/bachtis/ZMC.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>80&&Z_mass<120&&MuPosGenStatus1_pt>0&&MuNegGenStatus1_pt>0&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/ZDATA.root",1,0)
p.processTree("/tmp/bachtis/ZDATA.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>80&&Z_mass<120&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()


p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/YMC.root",0,0)
p.processTree("/tmp/bachtis/YMC.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>9.1&&Z_mass<9.7&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/YGEN.root",0,1)
p.processTree("/tmp/bachtis/YMC.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>9.1&&Z_mass<9.7&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5&&MuPosGenStatus1_pt>0&&MuNegGenStatus1_pt>0")
p.write()

p=ROOT.NtupleProcessorLuca("/data/bachtis/CALIB/7TeV/YDATA.root",1,0)
p.processTree("/tmp/bachtis/YDATA.root","MuPos_pt>0&&MuNeg_pt>0&&MuPos_charge+MuNeg_charge==0&&Z_mass>9.1&&Z_mass<9.7&&abs(MuPos_eta)<2.5&&abs(MuNeg_eta)<2.5")
p.write()



