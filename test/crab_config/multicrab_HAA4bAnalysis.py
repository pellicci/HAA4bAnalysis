from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'crab_projects/samples/'

config.section_('JobType')
config.JobType.psetName = 'run_HAA4bAnalysis.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['MC_Recent_25ns_2015.root','pileUpData_fromJson.root'] #data files for PileUp reweighting
config.JobType.outputFiles = ['HAA4bAnalysis_output.root']
config.JobType.pyCfgParams = ['runningOnData=False']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.section_('Site')
config.Site.storageSite = 'T2_IT_Legnaro'
#config.Site.storageSite = 'T2_IN_TIFR'

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

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'HAA4bAnalysis_ttbar'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_ttbarW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_ttbarZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_SingleTop_tW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_SingleAntiTop_tW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_WJetsToLNu'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_WJetsToQQ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WJetsToQQ_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_DY_10_50'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_DY_50'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT100to200'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT200to300_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join() 

    config.General.requestName = 'HAA4bAnalysis_QCD_HT200to300_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join() 

    config.General.requestName = 'HAA4bAnalysis_QCD_HT300to500_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT300to500_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT500to700_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT500to700_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT700to1000_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT700to1000_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT1000to1500_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT1000to1500_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT1500to2000_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = ' /QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT1500to2000_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = ' /QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT2000toInf_1'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT2000toInf_1_2'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_ZZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_WW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WWTo4Q_13TeV-powheg/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_WZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
