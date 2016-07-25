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

config.section_('Data')
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
#config.Data.outLFNDirBase = '/store/user/aashah/'
config.Data.publication = False

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

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################

    config.General.requestName = 'HAA4bAnalysis_ttbar'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12_ext4-v1/MINIAODSIM' ##check it for _ext4 which is an extra part 
   # config.Data.inputDataset = '/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_ttbarW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


    config.General.requestName = 'HAA4bAnalysis_ttbarZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()



#    config.General.requestName = 'HAA4bAnalysis_QCD_15_30'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()
#
#    config.General.requestName = 'HAA4bAnalysis_QCD_30_50'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()

#    config.General.requestName = 'HAA4bAnalysis_QCD_50_80'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#   p.join()

#   config.General.requestName = 'HAA4bAnalysis_QCD_80_120'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()

#   config.General.requestName = 'HAA4bAnalysis_QCD_120_170'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()

#   config.General.requestName = 'HAA4bAnalysis_QCD_170_300'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()
#
#   config.General.requestName = 'HAA4bAnalysis_QCD_300_470'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()
#
#   config.General.requestName = 'HAA4bAnalysis_QCD_470_600'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()

#   config.General.requestName = 'HAA4bAnalysis_QCD_600_800'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()

#   config.General.requestName = 'HAA4bAnalysis_QCD_800_1000'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()

#   config.General.requestName = 'HAA4bAnalysis_QCD_1000_1400'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()
#
#   config.General.requestName = 'HAA4bAnalysis_QCD_1400_1800'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()
#
#   config.General.requestName = 'HAA4bAnalysis_QCD_1800_2400'
#   config.Data.unitsPerJob = 5
#   config.Data.inputDataset = '/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#   p = Process(target=submit, args=(config,))
#   p.start()
#   p.join()

    config.General.requestName = 'HAA4bAnalysis_WJetsToLNu'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

#    config.General.requestName = 'HAA4bAnalysis_DY_5_50'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()

#    config.General.requestName = 'HAA4bAnalysis_DY_50'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()


    config.General.requestName = 'HAA4bAnalysis_QCD_HT100to200'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT200to300'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join() 

    config.General.requestName = 'HAA4bAnalysis_QCD_HT300to500'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT500to700'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT700to1000'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT1000to1500'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT1500to2000'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_HT2000toInf'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
#===========================================Another add

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P20to30'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-20to30_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P30to50'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-30to50_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P50to80'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-50to80_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P80to120'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-80to120_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P120to170'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-120to170_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P170to300'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P300to470'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P470to600'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P600to800'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-600to800_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P800to1000'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_QCD_MuEnriched_P1000toInf'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCUETP8M1_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

#======================================================
    config.General.requestName = 'HAA4bAnalysis_SingleTop_tW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()
    p.join()

    config.General.requestName = 'HAA4bAnalysis_SingleAntiTop_tW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

#    config.General.requestName = 'HAA4bAnalysis_ZZ'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/ZZ_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()

    config.General.requestName = 'HAA4bAnalysis_ZZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


#    config.General.requestName = 'HAA4bAnalysis_WW'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/WW_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()

    config.General.requestName = 'HAA4bAnalysis_WW'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WWTo4Q_13TeV-powheg/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()

#    config.General.requestName = 'HAA4bAnalysis_WZ'
#    config.Data.unitsPerJob = 5
#    config.Data.inputDataset = '/WZ_TuneCUETP8M1_13TeV-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
#    p = Process(target=submit, args=(config,))
#    p.start()
#    p.join()

    config.General.requestName = 'HAA4bAnalysis_WZ'
    config.Data.unitsPerJob = 5
    config.Data.inputDataset = '/WZJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    p = Process(target=submit, args=(config,))
    p.start()
    p.join()


