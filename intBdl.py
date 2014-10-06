#

import FWCore.ParameterSet.Config as cms

process = cms.Process("MAGNETICFIELDTEST")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

#process.load("MagneticField.Engine.uniformMagneticField_cfi")
#process.UniformMagneticFieldESProducer.ZFieldInTesla = 3.8

process.load("MagneticField.Engine.volumeBasedMagneticField_090322_2pi_scaled_cfi")
#process.load("MagneticField.ParametrizedEngine.parametrizedMagneticField_OAE_3_8T_cfi")
#process.load("MagneticField.ParametrizedEngine.parametrizedMagneticField_PolyFit2D_cfi")
#process.load("MagneticField.ParametrizedEngine.parametrizedMagneticField_PolyFit3D_cfi")
#process.load("MagneticField.Engine.volumeBasedMagneticField_120812_cfi")
#process.load("MagneticField.Engine.volumeBasedMagneticField_120812_largeYE4_cfi")

#process.VolumeBasedMagneticFieldESProducer.useParametrizedTrackerField = False
process.VolumeBasedMagneticFieldESProducer.scalingVolumes = cms.vint32()
process.VolumeBasedMagneticFieldESProducer.scalingFactors = cms.vdouble()

# process.MessageLogger = cms.Service("MessageLogger",
#     categories   = cms.untracked.vstring("MagneticField"),
#     destinations = cms.untracked.vstring("cout"),
#     cout = cms.untracked.PSet(  
#     noLineBreaks = cms.untracked.bool(True),
#     threshold = cms.untracked.string("INFO"),
#     INFO = cms.untracked.PSet(
#       limit = cms.untracked.int32(0)
#     ),
#     WARNING = cms.untracked.PSet(
#       limit = cms.untracked.int32(0)
#     ),
#     MagneticField = cms.untracked.PSet(
#      limit = cms.untracked.int32(10000000)
#     )
#   )
# )

process.integral  = cms.EDAnalyzer("intBdl")
process.p1 = cms.Path(process.integral)

