import FWCore.ParameterSet.Config as cms

HAA4bAnalysis = cms.EDAnalyzer('HAA4bAnalysis',
                               jets = cms.InputTag("slimmedJets"),
                               BTagAlgo = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                               minPt_high = cms.double(50.),
                               minPt_low = cms.double(20.),
                               minCSV3 = cms.double(0.605),
                               runningOnData = cms.bool(True),  #Must be TRUE when running on data
                               pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"), #New Stuff 
                               bsCollection = cms.InputTag("offlineBeamSpot"),
                               PileupSrc = cms.InputTag("slimmedAddPileupInfo")
)
