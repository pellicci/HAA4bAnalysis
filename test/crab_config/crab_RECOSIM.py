from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4b_Pythia8_RECOSIM_80XV1_MH800_MA300'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'HAA4b_13TeV_pythia8_RECO_cfg.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/HAA4b_GENSIM_80XV1_MH800_MA300/pellicci-HAA4b_GENSIM_80XV1_MH800_MA300-619a232e452025736da7dd1e9cea282f/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'HAA4b_RECOSIM_80XV1_MH800_MA300'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

