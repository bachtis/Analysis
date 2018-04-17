import FWCore.ParameterSet.Config as cms



process = cms.Process("Demo")



### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use#
    fileNames = cms.untracked.vstring(
'file:/tmp/bachtis/test.root'
    )
)

process.analysis = cms.EDAnalyzer('MuonCalibAnalyzer',
                                  muons = cms.InputTag("slimmedMuons"),
                                  met = cms.InputTag("slimmedMETs"),
                                  genParticles = cms.InputTag("packedGenParticles"),
                                  vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                  isOnia = cms.bool(True)
)


process.p = cms.Path(process.analysis)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')


#import PhysicsTools.PythonAnalysis.LumiList as LumiList
#jsonFile='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_27Jan2016ReReco_Collisions15_ZeroTesla_25ns_JSON_MuonPhys.txt'
#myLumis = LumiList.LumiList(filename = jsonFile).getCMSSWString().split(',')
#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange()
#process.source.lumisToProcess.extend(myLumis)


