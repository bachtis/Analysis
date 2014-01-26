import copy
import os 
import CMGTools.RootTools.fwlite.Config as cfg
from CMGTools.Production.datasetToSource import datasetToSource

hcalAnalyzer = cfg.Analyzer(
    'HCALShowerAnalyzer',
    )




tree = cfg.Analyzer(
    'HCALShowerTree'
)    

sequence=[
    hcalAnalyzer,
    tree
    ]


sequence = cfg.Sequence(sequence)


def getFiles(dataset, user, pattern):
    # print 'getting files for', dataset,user,pattern
    ds = datasetToSource( user, dataset, pattern, False )
    files = ds.fileNames
    return ['root://eoscms//eos/cms%s' % f for f in files]

upgrade = cfg.MCComponent(
    dataset='/Upgrade/SinglePi_1/GEN-SIM-DIGI-RECO',
    name = 'upgrade',
    files = getFiles('/Upgrade/SinglePi_20/GEN-SIM-DIGI/RECO','bachtis','.*root')+ \
            getFiles('/Upgrade/SinglePi_2/GEN-SIM-DIGI/RECO','bachtis','.*root')+ \
            getFiles('/Upgrade/SinglePi_5/GEN-SIM-DIGI/RECO','bachtis','.*root')+ \
            getFiles('/Upgrade/SinglePi_10/GEN-SIM-DIGI/RECO','bachtis','.*root')+ \
            getFiles('/Upgrade/SinglePi_20/GEN-SIM-DIGI/RECO','bachtis','.*root')+ \
            getFiles('/Upgrade/SinglePi_50/GEN-SIM-DIGI/RECO','bachtis','.*root') ,
#    files = ['file:reco.root'],
    xSection = 1,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1,
    isMC = True,
    isData = False,
    splitFactor =1
)



selectedComponents=[upgrade]
    
config = cfg.Config( components = selectedComponents,
                     sequence = sequence )
