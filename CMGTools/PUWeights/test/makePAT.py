## import skeleton process
from PhysicsTools.PatAlgos.patTemplate_cfg import *


## ------------------------------------------------------
#  NOTE: you can use a bunch of core tools of PAT to
#  taylor your PAT configuration; for a few examples
#  uncomment the lines below
## ------------------------------------------------------



from PhysicsTools.PatAlgos.tools.pfTools import *
usePFIso(process,"")

## remove certain objects from the default sequence
# removeAllPATObjectsBut(process, ['Muons'])
removeSpecificPATObjects(process, ['Electrons', 'Jets', 'Taus','Photons'])




process.skim = cms.EDFilter("MuMuSKIM",
                            process=cms.string("HLT")
)                            


## remove MC matching from the default sequence
from PhysicsTools.PatAlgos.tools.coreTools import *
removeMCMatching(process, ['Muons'])
runOnData(process)

process.patMuons.embedTrack=cms.bool(True)
process.patMuons.embedCaloMETMuonCorrs = cms.bool(False)
process.patMuons.embedTcMETMuonCorrs   = cms.bool(False)


#process.weightedNeutrals = cms.EDProducer('PFNeutralPUWeighter')
process.particleFlowWeighted = cms.EDProducer('PFPUMitigator')

## let it run
process.p = cms.Path(
    process.skim*
    process.patDefaultSequence*
    process.particleFlowWeighted
    )




## ------------------------------------------------------
#  In addition you usually want to change the following
#  parameters:
## ------------------------------------------------------
#
#   process.GlobalTag.globaltag =  ...    ##  (according to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions)
#                                         ##
from PhysicsTools.PatAlgos.patInputFiles_cff import filesRelValProdTTbarAODSIM
process.source.fileNames = filesRelValProdTTbarAODSIM


from CMGTools.Production.datasetToSource import *

datasetInfo = (
    'CMS',
    '/DoubleMu/Run2012A-22Jan2013-v1/AOD',
    '.*root')
process.source = datasetToSource(
    *datasetInfo
    )



#                                         ##
process.maxEvents.input = 10000
#                                         ##
process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_*')
#process.out.outputCommands.append('keep *_weightedNeutrals_*_*')
process.out.outputCommands.append('keep *_particleFlow_*_*')
process.out.outputCommands.append('keep *_particleFlowWeighted_*_*')
#                                         ##
process.out.fileName = 'patTuple_standard.root'
#                                         ##
#   process.options.wantSummary = False   ##  (to suppress the long output at the end of the job)

