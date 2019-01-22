import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras

process = cms.Process('HGCalSimTest',eras.Run2_2017,eras.run3_GEM)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.Geometry.GeometryExtended2023D23Reco_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('RecoLocalMuon.GEMRecHit.gemLocalReco_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v10', '')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
#process.maxEvents.input = cms.untracked.int32(10)
# Input source
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.skipEvents = cms.untracked.uint32(0)

process.source.fileNames.append('/store/user/yekang/me0/tenMu_default/tenMu_GEN-SIM-DIGI_060.root')
#from glob import glob
#process.source.fileNames.append('file:/cms/ldap_home/yckang/me0/CMSSW_10_4_0/src/MuonPerformance/MuonAnalyser/test/tenMu_GEN-SIM-DIGI.root')

#fname = 'singleMuon.txt'
#f = open(fname)
#for line in f:
#    process.source.fileNames.append(line)

process.options = cms.untracked.PSet(
)

process.TFileService = cms.Service("TFileService",fileName = cms.string("histo.root"))

process.HGCalSimTest = cms.EDAnalyzer('HGCalSimTest',
    me0Digis = cms.InputTag("simMuonME0Digis"),
    me0Segments = cms.InputTag("me0Segments"),
)
process.p = cms.Path(process.HGCalSimTest)

process.schedule = cms.Schedule(process.p)
