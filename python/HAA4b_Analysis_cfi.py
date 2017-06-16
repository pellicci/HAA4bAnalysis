import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
#####
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "import runningOnData here")
options.parseArguments()
########
if not options.runningOnData:
     jecLevels = ([
      'input_data/Summer16_23Sep2016V4_MC_L1FastJet_AK4PFchs.txt',   #Updated files, but still does have some problems
      'input_data/Summer16_23Sep2016V4_MC_L1RC_AK4PFchs.txt',
      'input_data/Summer16_23Sep2016V4_MC_L2L3Residual_AK4PFchs.txt',
      'input_data/Summer16_23Sep2016V4_MC_L2Relative_AK4PFchs.txt'
      'input_data/Summer16_23Sep2016V4_MC_L3Absolute_AK4PFchs.txt'
    ])
if options.runningOnData:
    jecLevels = ([
     'input_data/Summer16_23Sep2016BCDV4_DATA_L1FastJet_AK4PFchs.txt',  #For Samples B, C, and D
     'input_data/Summer16_23Sep2016BCDV4_DATA_L1RC_AK4PFchs.txt',
     'input_data/Summer16_23Sep2016BCDV4_DATA_L2L3Residual_AK4PFchs.txt',
     'input_data/Summer16_23Sep2016BCDV4_DATA_L2Relative_AK4PFchs.txt',
     'input_data/Summer16_23Sep2016BCDV4_DATA_L2Residual_AK4PFchs.txt',
     'input_data/Summer16_23Sep2016BCDV4_DATA_L3Absolute_AK4PFchs.txt'  #,
     #'input_data/Summer16_23Sep2016EFV4_DATA_L1FastJet_AK4PFchs.txt',               #For Samples E and F
     #'input_data/Summer16_23Sep2016EFV4_DATA_L1RC_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016EFV4_DATA_L2L3Residual_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016EFV4_DATA_L2Relative_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016EFV4_DATA_L2Residual_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016EFV4_DATA_L3Absolute_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016GV4_DATA_L1FastJet_AK4PFchs.txt',              #For Sample G
     #'input_data/Summer16_23Sep2016GV4_DATA_L1RC_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016GV4_DATA_L2L3Residual_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016GV4_DATA_L2Relative_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016GV4_DATA_L2Residual_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016GV4_DATA_L3Absolute_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016HV4_DATA_L1FastJet_AK4PFchs.txt',              #For Sample H
     #'input_data/Summer16_23Sep2016HV4_DATA_L1RC_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016HV4_DATA_L2L3Residual_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016HV4_DATA_L2Relative_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016HV4_DATA_L2Residual_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016HV4_DATA_L3Absolute_AK4PFchs.txt',
     #'input_data/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt'              
    ])
##
##-------------------- User analyzer  --------------------------------
HAA4bAnalysis = cms.EDAnalyzer('HAA4bAnalysis',
                               jets             = cms.InputTag("slimmedJets"),
                               globaljets       = cms.InputTag("slimmedJets"),
                               genParticles     = cms.InputTag("prunedGenParticles"),
                               met              = cms.InputTag("slimmedMETs"),
                               globalmet        = cms.InputTag("slimmedMETs"), #In case we apply a cut on the 'previous case'
                               BTagAlgo         = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                               minPt_high       = cms.double(50.),
                               minPt_low        = cms.double(20.),
                               minCSV           = cms.double(0.460),
                               runningOnData    = cms.bool(False),  #Controlled Externally depending on data or MC
                               pvCollection     = cms.InputTag("offlineSlimmedPrimaryVertices"), 
                               bsCollection     = cms.InputTag("offlineBeamSpot"),
                               PileupSrc        = cms.InputTag("slimmedAddPileupInfo"),
                               rhoSrc           = cms.InputTag("fixedGridRhoAll"),
			       jecPayloadNames  = cms.vstring(jecLevels),
 			       jecUncName       = cms.string("input_data/Summer16_23Sep2016V4_MC_Uncertainty_AK4PFchs.txt")   #MC
 			       #jecUncName       = cms.string("input_data/Summer16_23Sep2016BCDV4_DATA_Uncertainty_AK4PFchs.txt") #B,C, and D
 			       #jecUncName       = cms.string("input_data/Summer16_23Sep2016EFV4_DATA_Uncertainty_AK4PFchs.txt")  #E and F
 			       #jecUncName       = cms.string("input_data/Summer16_23Sep2016GV4_DATA_Uncertainty_AK4PFchs.txt")   #G
 			       #jecUncName       = cms.string("input_data/Summer16_23Sep2016HV4_DATA_Uncertainty_AK4PFchs.txt")   #H
)
