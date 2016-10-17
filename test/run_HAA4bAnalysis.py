import FWCore.ParameterSet.Config as cms
process = cms.Process("USER")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load('Configuration.StandardSequences.Services_cff') 
from Configuration.AlCa.GlobalTag import GlobalTag

#JEC
process.load('JetMETCorrections.Configuration.DefaultJEC_cff') 
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
from JetMETCorrections.Configuration.DefaultJEC_cff import *
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('runningOnData',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "PU config flag")
options.parseArguments()

#Input source
if options.runningOnData: 
   process.GlobalTag = GlobalTag(process.GlobalTag, '76X_dataRun2_v15') #which conditions to use
   print "Data Sample will be taken as input for check up of the code working "
   inputFiles="root://xrootd.unl.edu//store/data/Run2015D/BTagCSV/MINIAOD/16Dec2015-v1/50000/00AF8EB4-70AB-E511-9271-00266CFAE7AC.root"
else:
   process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')
   print "MC Sample will be taken as input for check up of the code working "
   inputFiles="root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv2/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/00000/0AAD5298-DBB8-E511-8527-003048D2BD8E.root" 

process.source = cms.Source ("PoolSource",
                             #fileNames = cms.untracked.vstring (inputFiles),
                             fileNames = cms.untracked.vstring('file:miniAOD-prod_new_PAT.root'), #When running on crab
)

# Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("HAA4bAnalysis_output.root")
)

#Put a loose selection on the b-jets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
process.selectedJets = cms.EDFilter("PATJetSelector",
                                     src = cms.InputTag("slimmedJets"),
                                     cut = cms.string("pt > 30 && abs(eta) < 2.5")
                                     #cut = cms.string('bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > 0.605')
                                     )

process.load("HiggsAnalysis.HAA4bAnalysis.HAA4b_Analysis_cfi")
process.HAA4bAnalysis.runningOnData = options.runningOnData

# Applying Jet Energy Corrections
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if not process.HAA4bAnalysis.runningOnData:      #This loop is for MC
        print "Jet Energy Corrections on Monte Carlo will be applied "
	updateJetCollection(
	   process,
	   jetSource = cms.InputTag('slimmedJets'),
	   labelName = 'UpdatedJEC',
	   jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None') # Do not forget 'L2L3Residual' on data!
	)
else:                                             #this loop is for data
        print "Jet Energy Corrections on Data will be applied "
	updateJetCollection(
	   process,
	   jetSource = cms.InputTag('slimmedJets'),
	   labelName = 'UpdatedJEC',
	   jetCorrections = ('AK5PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')# Added 'L2L3Residual'for data!
	)	
process.JetCorr = process.selectedJets.clone(src=cms.InputTag("updatedPatJetsUpdatedJEC"), cut = cms.string("pt > 30 && abs(eta) < 2.5"))

#process.HAA4bAnalysis.jets = cms.InputTag("selectedJets")         #Uncomment only when you need jets without energy corrections
#process.HAA4bAnalysis.jets = cms.InputTag("updatedPatJetsUpdatedJEC") #Pass on updated jets (jets with corrected energy) without any cut
process.HAA4bAnalysis.jets = cms.InputTag("JetCorr") #Pass on updated Jets with the required cuts.

process.HAA4bAnalysis.minPt_low = cms.double(30.)

import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.trigger_filter = hlt.triggerResultsFilter.clone()
process.trigger_filter.triggerConditions = cms.vstring('HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v*', 'HLT_QuadJet45_TripleBTagCSV0p67_v*')
process.trigger_filter.hltResults = cms.InputTag( "TriggerResults", "", "HLT" )
process.trigger_filter.l1tResults = cms.InputTag("")
process.trigger_filter.throw = cms.bool( False )

#process.seq = cms.Path(process.trigger_filter* process.selectedJets * process.HAA4bAnalysis ) #Only for Uncorrected Jets
process.seq = cms.Path(process.trigger_filter * process.selectedJets * process.patJetCorrFactorsUpdatedJEC*process.updatedPatJetsUpdatedJEC * process.JetCorr* process.HAA4bAnalysis )
process.schedule = cms.Schedule(process.seq)
