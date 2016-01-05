import FWCore.ParameterSet.Config as cms

HZZ4bAnalysis = cms.EDAnalyzer('HAA4bAnalysis',
                               jets = cms.untracked.InputTag("selectedPatJetsAK5PFCHS"),
                               BTagAlgo = cms.untracked.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                               minPt = cms.untracked.double(30.)
)
