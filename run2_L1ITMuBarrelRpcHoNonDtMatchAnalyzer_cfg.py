import FWCore.ParameterSet.Config as cms
process = cms.Process('L1ITMU')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2015_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1ITMuonBarrelRpcHoAnalyzerNonMatchedDT_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1ITMuTriggerPrimitiveProducer_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1CSCTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1DTTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.L1RPCTFTrackConverter_cfi')
process.load('L1Trigger.L1IntegratedMuonTrigger.MBLTProducer_cfi')

process.load('L1Trigger.L1IntegratedMuonTrigger.MBLTProducer_cfi')

from Configuration.AlCa.autoCond import autoCond
process.GlobalTag.globaltag = autoCond['run2_mc'] 

#infile=['/store/user/psaxena/L1Trigger/HOUpgrade/Generation/SingleMuonGun/SingleMuMinus_Fall14_FlatPt-0to200_MCRUN2_72_V3_GEN_SIM_DIGI_RECO_L1/150127_084421/0000/SingleMuMinus_Fall14_FlatPt-0to200_MCRUN2_72_V1_GEN_SIM_DIGI_RECO_L1_406.root']

#infile=['/store/user/psaxena/L1Trigger/HOUpgrade/Generation/SingleMuonGun/GEN-SIM-DIGI/SingleMuMinus_Winter15_FlatPt-0to200_MCRUN2_72_V3_GEN_SIM_DIGI/150316_134931/0000/SingleMuPlus_Fall14_FlatPt-1to140_PRE_LS172_V15_GEN_SIM_DIGI_69.root']

infile=['/store/user/psaxena/L1Trigger/HOUpgrade/Generation/SingleMuonGun/SingleMuMinus_Fall14_FlatPt-0to200_MCRUN2_72_V3_GEN_SIM_DIGI_RECO_L1/150127_084421/0002/SingleMuMinus_Fall14_FlatPt-0to200_MCRUN2_72_V1_GEN_SIM_DIGI_RECO_L1_2588.root']

process.source = cms.Source(
    'PoolSource',
    #    fileNames = cms.untracked.vstring(infile)
    #    fileNames = cms.untracked.vstring("file:SingleMuMinus_Winter15_FlatPt-0to200_MCRUN2_72_V1_GEN_SIM_DIGI_RECO_L1.root")
    fileNames = cms.untracked.vstring("file:SingleMuPlus_Fall14_FlatPt-3to200_PRE_LS172_V15_GEN_SIM_DIGI.root")
    #   fileNames = cms.untracked.vstring("/store/user/psaxena/L1Trigger/HOUpgrade/Generation/SingleMuonGun/GEN-SIM-DIGI/SingleMuMinus_Winter15_FlatPt-0to200_MCRUN2_72_V3_GEN_SIM_DIGI/150316_134931/0000/SingleMuPlus_Fall14_FlatPt-1to140_PRE_LS172_V15_GEN_SIM_DIGI_508.root")
    )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.L1MuBarrelRPCHOAnalyzer = process.L1ITMuBarrelRpcHoNonMatchDT.clone(
    minDistForRpcHOMatch = cms.double( 0.5 ), 
    )


process.L1ITMUSequence = cms.Sequence( process.L1ITMuTriggerPrimitives +
                                       process.L1CSCTFTrackConverter   +
                                       process.L1DTTFTrackConverter    +
                                       process.L1RPCTFTrackConverters  +
                                       process.MBLTProducer            +
                                       process.L1MuBarrelRPCHOAnalyzer)
#process.MBLTProducer)
#                                       process.L1ITMuBarrelRPCHOAnalyzer )
process.L1ITMUPath = cms.Path(process.L1ITMUSequence)

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("L1ITRpcHOalgo.root"),
      closeFileFast = cms.untracked.bool(True)
  )


#outCommands = cms.untracked.vstring('drop *')
#outCommands.append('keep *_genParticles_*_*')
#outCommands.append('keep *_simCsctfDigis_*_*')
#outCommands.append('keep *_simDttfDigis_*_*')
#outCommands.append('keep *_simRpcTriggerDigis_*_*')
#outCommands.append('keep *_simMuonRPCDigis_*_*')
##outCommands.append('keep HOData*_simHcalDigis_*_*')
#outCommands.append('keep *_simHcalDigis_*_*')
#outCommands.append('keep *_simDtTriggerPrimitiveDigis_*_*')
#outCommands.append('keep *_simCscTriggerPrimitiveDigis_*_*')
#outCommands.append('keep *_L1ITMuTriggerPrimitives_*_*')
#outCommands.append('keep *_*Converter_*_*')
#outCommands.append('keep *_*Matcher_*_*')
#

#process.FEVTDEBUGoutput = cms.OutputModule(
#    "PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
##    outputCommands = outCommands,
#    fileName = cms.untracked.string('L1ITMuBarrelRpcHo.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('')
#    )
#)
#
#
#process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
#)
#
#process.outPath = cms.EndPath(process.FEVTDEBUGoutput)

process.schedule = cms.Schedule(process.L1ITMUPath)

