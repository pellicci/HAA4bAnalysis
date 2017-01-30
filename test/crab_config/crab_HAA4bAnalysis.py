from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4bAnalysis_Signal_H800_A300'
config.General.workArea = 'crab_projects/samples'

config.section_('JobType')
config.JobType.psetName = 'run_HAA4bAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MC_Recent_25ns_2015.root','pileUpData_fromJson.root'] #files for PU reweighting
config.JobType.outputFiles = ['HAA4bAnalysis_output.root']
config.JobType.pyCfgParams = ['runningOnData=False'] #Set to True if MC misses PileUP info


config.section_('Data')
config.Data.inputDataset = '/HAA4b_GENSIM_80XV1_MH800_MA300/pellicci-HAA4b_MINIAODSIM_80XV3_MH800_MA300-0e6df83a66ae3f10341eebd4053e4881/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
#config.Site.storageSite = 'T2_IN_TIFR'
