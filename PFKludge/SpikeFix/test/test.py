import FWCore.ParameterSet.Config as cms

process = cms.Process("FIX")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.load("RecoJets.JetProducers.ak5PFJets_cfi")

process.ak5PFJets.src = cms.InputTag('fixedPF')

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/tmp/bachtis/upgrade2.root'
    )
)

process.fixedPF = cms.EDProducer('SpikeFix')

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('/tmp/bachtis/out.root')
)
  
process.p = cms.Path(process.fixedPF+process.ak5PFJets)

process.e = cms.EndPath(process.out)
