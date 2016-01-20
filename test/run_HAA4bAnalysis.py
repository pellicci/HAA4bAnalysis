import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

#Rerun jet reco for AK5 and new Btagging
from HiggsAnalysis.HAA4bAnalysis.reRunAK5Jets_Btag_cfg import *
reSetJet(process)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:miniAOD-prod_PAT.root'),
    secondaryFileNames = cms.untracked.vstring()
)

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HAA4bAnalysis_output.root")
)


#Put a loose selection on the b-jets
getattr(process,'selectedPatJetsAK5PFCHS').cut = cms.string('bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89')

process.load("HiggsAnalysis.HAA4bAnalysis.HAA4b_Analysis_cfi")
process.analysis = cms.Path(process.HZZ4bAnalysis)

process.schedule = cms.Schedule(process.rejet , process.analysis)
