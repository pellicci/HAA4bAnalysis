from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4b_Pythia8_MINIAODSIM_80XV1_MH800_MA300'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'HAA4b_13TeV_pythia8_MINIAOD_cfg.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/HAA4b_GENSIM_80XV1_MH800_MA300/pellicci-HAA4b_RECOSIM_80XV1_MH800_MA300-8cf88868cf066f834b46562c885350ae/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'HAA4b_MINIAODSIM_80XV1_MH800_MA300'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

