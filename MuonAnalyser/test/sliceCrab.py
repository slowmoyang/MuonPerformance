from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'GEM-SliceTestBkg-2018D-RAW-7036'
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'RAW2STA_SliceTest_Bkg.py'

#config.Data.inputDataset = '/EGamma/Run2018C-PromptReco-v1/AOD'
config.Data.inputDataset = '/ZeroBias/Run2018D-v1/RAW'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = True
# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'CRAB3_GEMSliceTestBkgAnalysis_GEM_2018D-7036'

# These values only make sense for processing data
#    Select input data based on a lumi mask
#config.Data.lumiMask = '2018D_321069.txt' #  'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
#    Select input data based on run-ranges
#config.Data.runRange = '321067-321069'  # v1
config.Data.runRange = '321908,321909'
#config.Data.runRange = '319347,319348,319349'  # v2

# Where the output files will be transmitted to
config.Site.storageSite = 'T3_KR_KISTI'
