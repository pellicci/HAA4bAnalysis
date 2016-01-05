from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'HAA4b_Pythia8_AODSIM'
config.General.workArea = 'crab_projects'

config.section_('JobType')
config.JobType.psetName = 'HAA4b_13TeV_pythia8_FULLSIM_cfg.py'
config.JobType.pluginName = 'PrivateMC'
config.JobType.outputFiles = ['HAA4b_pythia8_FULLSIM.root']
config.JobType.inputFiles = ['./L1TriggerConfig/L1GtConfigProducers/data/Luminosity/startup/L1Menu_Collisions2015_25nsStage1_v5_L1T_Scales_20141121.xml']

config.section_('Data')
config.Data.outputPrimaryDataset = 'HAA4b'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 100
NJOBS = 100
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = True
config.Data.outputDatasetTag = 'HAA4b'

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'

