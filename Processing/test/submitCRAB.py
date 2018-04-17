import os


samples={}

samples['jpsi']=[

'/Charmonium/Run2017B-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017C-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017D-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017E-31Mar2018-v1/MINIAOD',
'/Charmonium/Run2017F-31Mar2018-v1/MINIAOD',
'/JPsiToMuMu_Pt20to100-pythia8-gun/RunIIFall17MiniAOD-RECOSIMstep_94X_mc2017_realistic_v10-v1/MINIAODSIM'

]



samples['upsilon']=[
'/MuOnia/Run2017B-31Mar2018-v1/MINIAOD',
'/MuOnia/Run2017C-31Mar2018-v1/MINIAOD',
'/MuOnia/Run2017E-31Mar2018-v1/MINIAOD',
'/MuOnia/Run2017F-31Mar2018-v1/MINIAOD',
'/UpsilonMuMu_UpsilonPt6_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer17MiniAOD-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/MINIAODSIM'

]

samples['z']=[
'/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD',
'/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD',
'/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAOD-94X_mc2017_realistic_v10-v1/MINIAODSIM'

]
/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt
json='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'


for resonance,slist in samples.iteritems():
    for i,s in enumerate(slist):
        filename="crabSubmit_"+resonance+"_"+str(i)+'.py'
        f=open(filename,'w')
        cfg="""
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()


config.General.requestName = '{name}'
config.General.workArea = 'CRAB'
config.General.transferOutputs = True
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '{cfgFile}'
config.JobType.outputFiles = ['muonTree.root']
config.Data.inputDataset = '{dataset}'
config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.lumiMask = '{json}'
config.Data.outLFNDirBase = '/store/user/bachtis/'
config.Data.publication = False
config.Site.storageSite = 'T3_US_FNALLPC'

""".format( name=resonance+"_"+str(i),cfgFile='runZ.py' if resonance=='z' else 'runOnia.py',dataset=s,json=json)

        f.write(cfg)
        f.close()



        os.system('crab submit --config='+filename)
