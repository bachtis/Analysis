

import ROOT 




p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JDATA.root",1,0)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/Charmonium_Run2015D_16Dec2015_25ns/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>5&&Onia_l1_pt<100&&Onia_l2_pt>5&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0")
p.write()


p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/YDATA.root",1,0)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/MuOnia_Run2015D_16Dec2015_25ns/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>9.1&&Onia_mll<9.7&&Onia_l1_pt>5&&Onia_l1_pt<100&&Onia_l2_pt>5&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0")
p.write()




p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JMC.root",0,0)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/JpsiToMuMuPt8/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>5&&Onia_l1_pt<100&&Onia_l2_pt>5&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()

p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JGEN.root",0,1)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/JpsiToMuMuPt8/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>5&&Onia_l1_pt<100&&Onia_l2_pt>5&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()



p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JMCLOW.root",0,0)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/pt3/JpsiToMuMuPt8/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>3&&Onia_l1_pt<100&&Onia_l2_pt>3&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()


p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JGENLOW.root",0,1)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/pt3/JpsiToMuMuPt8/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>3&&Onia_l1_pt<100&&Onia_l2_pt>3&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()






p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/YMC.root",0,0)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/UpsToMuMuPt6/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>9.1&&Onia_mll<9.7&&Onia_l1_pt>5&&Onia_l1_pt<100&&Onia_l2_pt>5&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()

p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/YGEN.root",0,1)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/UpsToMuMuPt6/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>9.1&&Onia_mll<9.7&&Onia_l1_pt>5&&Onia_l1_pt<100&&Onia_l2_pt>5&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()



p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/ZDATA.root",1,0,"z")
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Z/DoubleMuon_Run2015D_16Dec2015_25ns/twoLeptonTreeProducer/tree.root","(z_l1_charge+z_l2_charge)==0&&z_mll>70&&z_mll<120&&z_l1_pt>20&&z_l1_pt<200&&z_l2_pt>20&&z_l2_pt<200&&abs(z_l1_eta)<2.5&&abs(z_l2_eta)<2.5&&z_l1_ptErr>0&&z_l2_ptErr>0&&abs(z_l1_pdgId)==13")
p.write()


p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/ZMC.root",0,0,"z")
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Z/DYJetsToLL_LO_M50/twoLeptonTreeProducer/tree.root","(z_l1_charge+z_l2_charge)==0&&z_mll>70&&z_mll<120&&z_l1_pt>20&&z_l1_pt<200&&z_l2_pt>20&&z_l2_pt<200&&abs(z_l1_eta)<2.5&&abs(z_l2_eta)<2.5&&z_l1_ptErr>0&&z_l2_ptErr>0&&abs(z_l1_pdgId)==13&&z_l1_mcPt>0&&z_l2_mcPt>0")
p.write()


p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/ZGEN.root",0,1,"z")
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Z/DYJetsToLL_LO_M50/twoLeptonTreeProducer/tree.root","(z_l1_charge+z_l2_charge)==0&&z_mll>70&&z_mll<120&&z_l1_pt>20&&z_l1_pt<200&&z_l2_pt>20&&z_l2_pt<200&&abs(z_l1_eta)<2.5&&abs(z_l2_eta)<2.5&&z_l1_ptErr>0&&z_l2_ptErr>0&&abs(z_l1_pdgId)==13&&z_l1_mcPt>0&&z_l2_mcPt>0")
p.write()





