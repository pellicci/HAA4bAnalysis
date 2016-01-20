from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4bAnalysis_Signal_H500_A200'
config.General.workArea = 'crab_projects/samples'

config.section_('JobType')
config.JobType.psetName = 'run_HAA4bAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['HAA4bAnalysis_output.root']

config.section_('Data')
config.Data.inputDataset = '/HAA4b_AODSIM/pellicci-HAA4b_MINIAODSIM-fd5dd32fcb4957ded3ffe3767ed84cb2/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

