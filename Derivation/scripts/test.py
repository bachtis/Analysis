

import ROOT 









p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JMCLOW.root",0,0)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/pt3/JpsiToMuMuPt8/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>3&&Onia_l1_pt<100&&Onia_l2_pt>3&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0")
p.write()


p=ROOT.NtupleProcessor("/data/bachtis/CALIB/13TeV_76X/JGENLOW.root",0,1)
p.processTree("/afs/cern.ch/user/g/gpetrucc/w/TREES_HZZ4L_220116/2L/Onia/pt3/JpsiToMuMuPt8/twoLeptonTreeProducerOnia/tree.root","(Onia_l1_charge+Onia_l2_charge)==0&&Onia_mll>2.9&&Onia_mll<3.3&&Onia_l1_pt>3&&Onia_l1_pt<100&&Onia_l2_pt>3&&Onia_l2_pt<100&&abs(Onia_l1_eta)<2.5&&abs(Onia_l2_eta)<2.5&&Onia_l1_ptErr>0&&Onia_l2_ptErr>0&&Onia_l1_mcPt>0&&Onia_l2_mcPt>0")
p.write()









