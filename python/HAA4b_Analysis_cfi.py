import FWCore.ParameterSet.Config as cms

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "import runningOnData here")
options.parseArguments()

if not options.runningOnData:
     jecLevels = ([
      'PHYS14_25_V1_L1FastJet_AK4PFchs.txt',   #We might have to usee new updated text files in near future
      'PHYS14_25_V1_L2Relative_AK4PFchs.txt',
      'PHYS14_25_V1_L3Absolute_AK4PFchs.txt'
    ])
else :
    jecLevels = ([
     'PHYS14_25_V1_L1FastJet_AK4PFchs.txt',
     'PHYS14_25_V1_L2Relative_AK4PFchs.txt',
     'PHYS14_25_V1_L3Absolute_AK4PFchs.txt',
     'PHYS14_25_V1_L2L3Residual_AK4PFchs.txt'
     ])

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
                               #jecPayloadNames  = cms.vstring("PHYS14_25_V1_L3Absolute_AK4PFchs.txt","PHYS14_25_V1_L2Relative_AK4PFchs.txt" ), 
 			       jecUncName       = cms.string("PHYS14_25_V1_Uncertainty_AK4PFchs.txt")
)
