# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions 101X_dataRun2_Prompt_v10 -s RAW2DIGI,L1Reco,RECO --runUnscheduled --process reRECO --data --era Run2_2017 --eventcontent RECO --hltProcess reHLT --scenario pp --datatier RECO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2017 -n -1 --filein /store/data/Run2017G/SingleMuon/RAW/v1/000/306/801/00001/7ADC8664-3DCE-E711-ADEF-02163E019BF4.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('reRECO',eras.Run2_2018,eras.run3_GEM)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( #'file:/xrootd/store/data/Run2017G/SingleMuon/RAW/v1/000/306/826/00000/FE2B6447-A4CE-E711-8DD7-02163E019CD7.root'),
        # "file:/xrootd/store/data/Run2018C/SingleMuon/RAW/v1/000/319/347/00000/041ED888-4683-E811-B79A-FA163EF19885.root"
        "/store/user/iawatson/GemSimulation/step2/step2_000.root"
    ),
    secondaryFileNames = cms.untracked.vstring()
)

#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = 'cert_306824-306826_5TeV.txt').getVLuminosityBlockRange()
#print process.source.lumisToProcess

process.options = cms.untracked.PSet()
# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v7', '')  # Run on 2018 Data
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2018_realistic', '')  # Run on 2018 Simulation

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
# process.reconstruction_step = cms.Path(process.reconstruction)

# Seed generator
from RecoMuon.MuonSeedGenerator.standAloneMuonSeeds_cff import *

# Stand alone muon track producer
from RecoMuon.StandAloneMuonProducer.standAloneMuons_cff import *

# Beam Spot 
from RecoVertex.BeamSpotProducer.BeamSpot_cff import *

process.dump=cms.EDAnalyzer('EventContentAnalyzer')

process.reconstruction_step = cms.Path(
    (process.dt1DRecHits + process.dt4DSegments + process.dt1DCosmicRecHits + process.dt4DCosmicSegments + process.csc2DRecHits + process.cscSegments + process.rpcRecHits + process.gemRecHits + process.gemSegments) *
    offlineBeamSpot*standAloneMuonSeeds*process.standAloneMuons  # *process.dump
)

process.endjob_step = cms.EndPath(process.endOfProcess)

process.TFileService = cms.Service("TFileService",fileName = cms.string("gem_raw2aod.root"))

process.SliceTestAnalysis = cms.EDAnalyzer(
    'STASliceTestAnalysis',
    process.MuonServiceProxy,
    gemRecHits = cms.InputTag("gemRecHits"),
    muons = cms.InputTag("standAloneMuons"),
    vertexCollection = cms.InputTag("offlinePrimaryVertices"),
    lumiScalers = cms.InputTag("scalersRawToDigi"),
    gemDigis = cms.InputTag("muonGEMDigis", "AMCdata", "reRECO"), 
    amc13Event = cms.InputTag("muonGEMDigis", "AMC13Event", "reRECO"),
    amcData = cms.InputTag("muonGEMDigis", "AMCdata", "reRECO"),
    gebStatusCol = cms.InputTag("muonGEMDigis", "gebStatus", "reRECO"),
    vfatStatusCol = cms.InputTag("muonGEMDigis", "vfatStatus", "reRECO"), 
)
process.sliceTest = cms.Path(process.SliceTestAnalysis)

process.muonGEMDigis.unPackStatusDigis = cms.bool(True)
# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,
                                process.L1Reco_step,
                                process.reconstruction_step,
                                process.sliceTest,
                                process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2017 

#call to customisation function customisePostEra_Run2_2017 imported from Configuration.DataProcessing.RecoTLR
#process = customisePostEra_Run2_2017(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)

# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
# from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
# process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
