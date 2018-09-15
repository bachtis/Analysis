import FWCore.ParameterSet.Config as cms



process = cms.Process("Demo")



### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use#
    fileNames = cms.untracked.vstring(
        '/store/mc/RunIIFall15DR76/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/08DC4220-16A7-E511-AF59-1CC1DE19286E.root'
    )
)

process.analysis = cms.EDAnalyzer('MuonCalibAnalyzer',
                                  muons = cms.InputTag("muons"),
                                  met = cms.InputTag("pfMet"),
                                  genParticles = cms.InputTag("genParticles"),
                                  vertices = cms.InputTag("offlinePrimaryVertices"),
                                  isOnia = cms.bool(False)
)


process.p = cms.Path(process.analysis)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#jsonFile='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_27Jan2016ReReco_Collisions15_ZeroTesla_25ns_JSON_MuonPhys.txt'
#myLumis = LumiList.LumiList(filename = jsonFile).getCMSSWString().split(',')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)


