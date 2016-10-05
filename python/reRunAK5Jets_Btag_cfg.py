import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

def reSetJet(process):
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
#################################################################
    slimmedAddPileupInfo = cms.EDProducer(
       'PileupSummaryInfoSlimmer',
        #src = cms.InputTag('addPileupInfo'),
        keepDetailedInfoFor = cms.vint32(0)
    )
#######################################################
    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection

    ## Add PAT jet collection based on the above-defined ak5PFJetsCHS
    addJetCollection(
        process,
        labelName = 'AK5PFCHS',
        jetSource = cms.InputTag('ak5PFJetsCHS'),
        pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
        bsSource = cms.InputTag('offlineBeamSpot'),
        pfCandidates = cms.InputTag('packedPFCandidates'),
        svSource = cms.InputTag('slimmedSecondaryVertices'),
        btagDiscriminators = bTagDiscriminators,
        #jetCorrections = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),  #commented before adding jet corrections
        genJetCollection = cms.InputTag('ak5GenJetsNoNu'),
        genParticles = cms.InputTag('prunedGenParticles'),
        algo = 'AK',
        rParam = 0.5
    )

    getattr(process,'selectedPatJetsAK5PFCHS').cut = cms.string('pt > 20.')

    process.rejet = cms.Path(process.packedGenParticlesForJetsNoNu *
                             process.ak5GenJetsNoNu *
                             process.pfCHS *
                             process.ak5PFJetsCHS *
                             process.patJetCorrFactorsAK5PFCHS *
                             process.patJetPartons *
                             process.patJetFlavourAssociationAK5PFCHS *
                             process.patJetPartonMatchAK5PFCHS *
                             process.patJetGenJetMatchAK5PFCHS *
                             process.pfImpactParameterTagInfosAK5PFCHS *
                             process.pfInclusiveSecondaryVertexFinderTagInfosAK5PFCHS *
                             process.pfCombinedInclusiveSecondaryVertexV2BJetTagsAK5PFCHS *
                             process.patJetsAK5PFCHS *
                             process.selectedPatJetsAK5PFCHS)

    from PhysicsTools.PatAlgos.tools.pfTools import adaptPVs
    ## Adapt primary vertex collection and BeamSpot
    adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))
    adaptBSs(process, bsCollection=cms.InputTag('offlineBeamSpot'))
 
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
        process,
        jetSource = cms.InputTag('slimmedJets'),
        labelName = 'UpdatedJEC',
        jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None')  # Do not forget 'L2L3Residual' on data
     )

    process.load('Configuration.StandardSequences.Services_cff') #Added info from now onwards
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")

    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    process.jer = cms.ESSource("PoolDBESSource",
            CondDBSetup,
            toGet = cms.VPSet(
                # Resolution
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_Fall15_25nsV2_MC_PtResolution_AK5PFchs'),
                    label  = cms.untracked.string('AK5PFchs_pt')
                    ),

                # Scale factors
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_Fall15_25nsV2_MC_SF_AK4PFchs'),
                    label  = cms.untracked.string('AK5PFchs')
                    ),
                ),
            connect = cms.string('sqlite:Fall15_25nsV2_MC.db')
            )
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


