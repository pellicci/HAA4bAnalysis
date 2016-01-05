import FWCore.ParameterSet.Config as cms

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

#from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
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

process.p = cms.Path(
                     #process.packedGenParticlesForJetsNoNu *
                     #process.ak5GenJetsNoNu *
                     #process.pfCHS *
                     #process.ak5PFJetsCHS *
                     # process.patJetCorrFactorsAK5PFCHS *
                     # process.patJetPartons *
                     # process.patJetFlavourAssociationAK5PFCHS *
                     # process.patJetPartonMatchAK5PFCHS *
                     # process.patJetGenJetMatchAK5PFCHS *
                     # process.pfCombinedInclusiveSecondaryVertexV2BJetTagsAK5PFCHS *
                     # process.patJetsAK5PFCHS *
                     process.selectedPatJetsAK5PFCHS)

from PhysicsTools.PatAlgos.tools.pfTools import adaptPVs
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))
