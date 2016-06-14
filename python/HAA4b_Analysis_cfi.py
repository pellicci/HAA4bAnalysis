import FWCore.ParameterSet.Config as cms

HZZ4bAnalysis = cms.EDAnalyzer('HAA4bAnalysis',
                               jets = cms.InputTag("slimmedJets"),
                               BTagAlgo = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                               minPt1 = cms.double(50.),
                               minPt4 = cms.double(20.),
                               runningOnData=cms.bool(True),  #Must be TRUE when running on data
                               pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"), #New Stuff 
                               bsCollection = cms.InputTag("offlineBeamSpot"),
                               PileupSrc = cms.InputTag("slimmedAddPileupInfo"),
)
