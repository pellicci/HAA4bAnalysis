from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4b_Pythia8_MINIAODSIM'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'miniAOD-prod_PAT.py'
config.JobType.pluginName = 'Analysis'

config.section_('Data')
config.Data.inputDataset = '/HAA4b_AODSIM/pellicci-HAA4b_AODSIM-b5eca7e9a51ba8c3ff68efa7a1e9c4f1/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'HAA4b_MINIAODSIM'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

