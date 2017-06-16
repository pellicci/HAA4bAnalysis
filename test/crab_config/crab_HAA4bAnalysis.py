from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4bAnalysis_Signal_H800_A300'
config.General.workArea = 'crab_projects/samples'

config.section_('JobType')
config.JobType.psetName = 'run_HAA4bAnalysis.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.inputFiles = ['MCpileUp_25ns_Recent2016.root','pileUpHistogramFromjson_Nominal.root', 'pileUpHistogramFromjson_ScaleUp.root', 'pileUpHistogramFromjson_ScaleDown.root', 'CSVv2_Moriond17_B_H.csv','PHYS14_25_V1_L1FastJet_AK4PFchs.txt', 'PHYS14_25_V1_L2L3Residual_AK4PFchs.txt','PHYS14_25_V1_L2Relative_AK4PFchs.txt','PHYS14_25_V1_L3Absolute_AK4PFchs.txt', 'PHYS14_25_V1_Uncertainty_AK4PFchs.txt'] #data files for PileUp reweighting
config.JobType.inputFiles = ['MCpileUp_25ns_Recent2016.root','pileUpHistogramFromjson_Nominal.root', 'pileUpHistogramFromjson_ScaleUp.root', 'pileUpHistogramFromjson_ScaleDown.root', 'CSVv2_Moriond17_B_H.csv', 'input_data/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt', 'input_data/Summer16_23Sep2016V4_MC_L1RC_AK4PFchs.txt', 'input_data/Summer16_23Sep2016V4_MC_L2L3Residual_AK4PFchs.txt', 'input_data/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt', 'input_data/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt', 'input_data/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt'] #data files for PileUp reweighting
config.JobType.outputFiles = ['HAA4bAnalysis_output.root']
config.JobType.pyCfgParams = ['runningOnData=False'] #Set to True if MC misses PileUP info


config.section_('Data')
config.Data.inputDataset = '/HAA4b_GENSIM_80XV1_MH800_MA300/pellicci-HAA4b_MINIAODSIM_80XV3_MH800_MA300-0e6df83a66ae3f10341eebd4053e4881/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
#config.Site.storageSite = 'T2_IN_TIFR'
