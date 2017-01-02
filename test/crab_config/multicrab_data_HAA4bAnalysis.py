from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/data/'

config.section_('JobType')
config.JobType.psetName = '../run_HAA4bAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['HAA4bAnalysis_output.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.pyCfgParams = ['runningOnData=True']

config.section_('Data')
#config.Data.lumiMask = 'json/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_v2.txt'
config.Data.lumiMask = '../json/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.unitsPerJob = 50
config.Data.inputDBS = 'global' #'http://cmsdbsprod.cern.ch/cms_dbs_prod_global/servlet/DBSServlet'
config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/aashah/'
config.section_('User')
config.section_('Site')
#config.Site.storageSite = 'T2_IT_Legnaro'
config.Site.storageSite = 'T2_IN_TIFR'

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException
    from multiprocessing import Process

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #########    From now on that's what users should modify: this is the a-la-CRAB2 configuration part.
        
 
#    config.General.requestName = 'HAA4bAnalysis_BTagCSV_B2'
#    config.Data.inputDataset = '/BTagCSV/Run2016B-23Sep2016-v2/MINIAOD'
#    config.Data.unitsPerJob = 5
#    from multiprocessing import Process
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()

    config.General.requestName = 'HAA4bAnalysis_BTagCSV_B'
    config.Data.inputDataset = '/BTagCSV/Run2016B-23Sep2016-v3/MINIAOD'
    config.Data.unitsPerJob = 5
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


    config.General.requestName = 'HAA4bAnalysis_BTagCSV_C'
    config.Data.inputDataset = '/BTagCSV/Run2016C-23Sep2016-v1/MINIAOD'
    config.Data.unitsPerJob = 5
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()



    config.General.requestName = 'HAA4bAnalysis_BTagCSV_D'
    config.Data.inputDataset = '/BTagCSV/Run2016D-23Sep2016-v1/MINIAOD'
    config.Data.unitsPerJob = 5
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_BTagCSV_E'
    config.Data.inputDataset = '/BTagCSV/Run2016E-23Sep2016-v1/MINIAOD'
    config.Data.unitsPerJob = 5
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_BTagCSV_F'
    config.Data.inputDataset = '/BTagCSV/Run2016F-23Sep2016-v1/MINIAOD'
    config.Data.unitsPerJob = 5
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_BTagCSV_G'
    config.Data.inputDataset = '/BTagCSV/Run2016G-23Sep2016-v1/MINIAOD'
    config.Data.unitsPerJob = 5
    from multiprocessing import Process
    p = Process(target=submit, args=(config,))
    p.start()
