import FWCore.ParameterSet.Config as cms

HZZ4bAnalysis = cms.EDAnalyzer('HAA4bAnalysis',
                               jets = cms.InputTag("selectedPatJetsAK5PFCHS"),
                               BTagAlgo = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
                               minPt1 = cms.double(50.),
                               minPt4 = cms.double(20.)
)
