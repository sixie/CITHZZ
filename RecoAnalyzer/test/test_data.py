# Auto generated configuration file
# using: 
# Revision: 1.149 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: promptCollisionReco -s RAW2DIGI,L1Reco,RECO,DQM,ALCA:SiStripCalZeroBias --datatier RECO --eventcontent RECO --conditions CRAFT09_R_V4::All --scenario pp --no_exec --data --magField AutoFromDBCurrent -n 100
import FWCore.ParameterSet.Config as cms

process = cms.Process('RecoAnalysis')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # print one event every

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/L1Reco_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/mc/Fall10/QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6/GEN-SIM-RECO/START38_V12-v1/0005/D85FC8E6-8DCC-DF11-9479-0024E8769B2C.root'
    'file:E0FB3164-CC52-E011-B693-E41F13181884.spring11dyeem20tunez2pythia6.root'
    #'/store/data/Run2011A/Photon/RECO/May10ReReco-v1/0005/029FBE15-1D7C-E011-BC71-001A64789E44.root'

    ) 

)


process.GlobalTag.globaltag = 'FT_R_42_V21B::All'


process.hltTrigReport = cms.EDAnalyzer( 'HLTrigReport',
     HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' )
)


from RecoJets.Configuration.RecoPFJets_cff import *
process.kt6PFJetsForRhoComputationEtaMaxDefault = kt6PFJets.clone(doRhoFastjet = True)

process.kt6PFJetsForRhoComputationEtaMax25 = kt6PFJets.clone(Rho_EtaMax = 2.5,doRhoFastjet = True)



process.recoAnalyzer = cms.EDAnalyzer("RecoAnalyzer",
                                      outputFile = cms.string('analysis.PhotonRun2011A-30Nov2011-v1AOD.160329-163869.root'),
                                      trigResultTag = cms.untracked.InputTag('TriggerResults','','HLT'),
                                           trigSummaryTag =  cms.untracked.InputTag('hltTriggerSummaryAOD','','HLT'),
                                          rhoCorrection    = cms.untracked.InputTag('kt6PFJetsForRhoComputationEtaMax25','rho'),
                                          rhoCorrection2    = cms.untracked.InputTag('kt6PFJetsForRhoComputationEtaMaxDefault','rho'),
                                          debugLevel = cms.int32(0)
                                      )

process.out_step = cms.EndPath(process.recoAnalyzer )

process.rhocomputeDefaultEtaMax = cms.Sequence(process.kt6PFJetsForRhoComputationEtaMax25*process.kt6PFJetsForRhoComputationEtaMaxDefault)
process.rhopath2 = cms.Path(process.rhocomputeDefaultEtaMax)

process.schedule = cms.Schedule(process.rhopath2,process.out_step)

