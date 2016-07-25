from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4bAnalysis_Signal_H800_A300'
config.General.workArea = 'crab_projects/samples'

config.section_('JobType')
config.JobType.psetName = 'run_HAA4bAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['HAA4bAnalysis_output.root']
config.JobType.pyCfgParams = ['runningOnData=False']

config.section_('Data')
#config.Data.inputDataset = '/HAA4b_AODSIM_76X_MH800_MA300/pellicci-HAA4b_MINIAODSIM_76X_MH800_MA300-b640eb3109575ebf90f337afac3d4f41/USER'
config.Data.inputDataset = '/HAA4b_AODSIM_76X_MH800_MA300/pellicci-HAA4b_MINIAODSIM_76X_MH800_MA300-b640eb3109575ebf90f337afac3d4f41/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/aashah/'
config.Data.publication = False

config.section_('Site')
#config.Site.storageSite = 'T2_IT_Legnaro'
config.Site.storageSite = 'T2_IN_TIFR'


