import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("USER")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

####################################################################################
#Rerun jet reco for AK5 and new Btagging
from FWCore.ParameterSet.VarParsing import VarParsing

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
## Define GenJets
from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
process.ak5GenJetsNoNu = ak5GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

## Select charged hadron subtracted packed PF candidates
process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
from RecoJets.JetProducers.ak5PFJets_cfi import ak5PFJets
## Define PFJetsCHS
process.ak5PFJetsCHS = ak5PFJets.clone(src = 'pfCHS', doAreaFastjet = True)

#################################################
## Remake PAT jets
#################################################

## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]

from PhysicsTools.PatAlgos.tools.jetTools import *

## Add PAT jet collection based on the above-defined ak5PFJetsCHS
addJetCollection(
    process,
    labelName = 'AK5PFCHS',
    jetSource = cms.InputTag('ak5PFJetsCHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 0.5
)

getattr(process,'selectedPatJetsAK5PFCHS').cut = cms.string('pt > 10')

from PhysicsTools.PatAlgos.tools.pfTools import adaptPVs
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))
####################################################################################

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
#getattr(process,'selectedPatJetsAK5PFCHS').cut = cms.string('bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.605')

from HiggsAnalysis.HAA4bAnalysis.HAA4b_Analysis_cfi import *

