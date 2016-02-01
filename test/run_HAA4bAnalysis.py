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
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:miniAOD-prod_new_PAT.root'),
    secondaryFileNames = cms.untracked.vstring()
)

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HAA4bAnalysis_output.root")
)


#Put a loose selection on the b-jets
getattr(process,'selectedPatJetsAK5PFCHS').cut = cms.string('bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.89')

process.load("HiggsAnalysis.HAA4bAnalysis.HAA4b_Analysis_cfi")
#process.analysis = cms.Path(process.HZZ4bAnalysis)

import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.trigger_filter = hlt.triggerResultsFilter.clone()
process.trigger_filter.triggerConditions = cms.vstring('HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v*', 'HLT_QuadJet45_TripleBTagCSV0p67_v*')
process.trigger_filter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
process.trigger_filter.l1tResults = cms.InputTag("")
process.trigger_filter.throw = cms.bool( False )
#process.FilterPath = cms.Path(process.trigger_filter)

#process.schedule = cms.Sequence(process.rejet * process.analysis )
process.seq = cms.Path(process.trigger_filter * process.HZZ4bAnalysis )

process.schedule = cms.Schedule(process.rejet , process.seq)
