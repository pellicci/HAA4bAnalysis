import FWCore.ParameterSet.Config as cms

HAA4bAnalysis = cms.EDAnalyzer('HAA4bAnalysis',
                               jets = cms.InputTag("slimmedJets"),
                               globaljets = cms.InputTag("slimmedJets"),
                               genParticles = cms.InputTag("prunedGenParticles"),
                               met = cms.InputTag("slimmedMETs"),
                               globalmet = cms.InputTag("slimmedMETs"), #In case we apply a cut on the 'previous case'
                               eventInfo = cms.untracked.bool(True),
                               BTagAlgo = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                               minPt_high = cms.double(50.),
                               minPt_low = cms.double(20.),
                               minCSV = cms.double(0.460),
                               runningOnData = cms.bool(False),  #Must be TRUE when running on data if running locally to test the code
                               pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"), #New Stuff 
                               bsCollection = cms.InputTag("offlineBeamSpot"),
                               PileupSrc = cms.InputTag("slimmedAddPileupInfo") # ,
)
