// -*- C++ -*-
//
// Package:    Temp/HAA4bAnalysis
// Class:      HAA4bAnalysis
// 
/**\class HAA4bAnalysis HAA4bAnalysis.cc Temp/HAA4bAnalysis/plugins/HAA4bAnalysis.cc

 Description: Selection of events for H -> AA -> 4b

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mario Pelliccioni
//         Created:  Tue, 05 Jan 2016 10:25:36 GMT


//system include files
#include <memory>
#include <algorithm>

//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <string>
#include <vector>
#include <fstream>

//CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

//Generator Level Information
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//Kinematic fitter
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

//Vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//PileUp inclusions
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//Jet Energy corrections and Uncertainties
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetCorrectionESSource.h"
#include "JetMETCorrections/Modules/interface/JetCorrectionESChain.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrectorCalculator.h"

//Jet Energy Resolutions and Uncertainties
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "CondFormats/DataRecord/interface/JetResolutionRcd.h"
#include "CondFormats/DataRecord/interface/JetResolutionScaleFactorRcd.h"
#include <boost/shared_ptr.hpp>
#include "FWCore/Framework/interface/ESHandle.h"
#include <boost/shared_ptr.hpp>

//BTag Scale factors from database
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

typedef math::XYZTLorentzVector LorentzVector;

#include "HAA4bAnalysis.h"

// constructors and destructor
HAA4bAnalysis::HAA4bAnalysis(const edm::ParameterSet& iConfig) :
  jets_            (iConfig.getParameter<edm::InputTag>("jets")),
  globaljets_      (iConfig.getParameter<edm::InputTag>("globaljets")), 
  genParticles_    (iConfig.getParameter<edm::InputTag>("genParticles")),
  met_             (iConfig.getParameter<edm::InputTag>("met")), //Met
  globalmet_       (iConfig.getParameter<edm::InputTag>("globalmet")),
  bdiscr_          (iConfig.getParameter<std::string>("BTagAlgo")),
  minPt_high_      (iConfig.getParameter<double>("minPt_high")),
  minPt_low_       (iConfig.getParameter<double>("minPt_low")),
  minCSV_          (iConfig.getParameter<double>("minCSV")),
  runningOnData_   (iConfig.getParameter<bool>("runningOnData")),   
  pvCollection_    (iConfig.getParameter<edm::InputTag>("pvCollection")),   
  bsCollection_    (iConfig.getParameter<edm::InputTag>("bsCollection")),  
  PileupSrc_       (iConfig.getParameter<edm::InputTag>("PileupSrc")),    
  rhoSrc_          (iConfig.getParameter<edm::InputTag>("rhoSrc")),
  jecPayloadNames_ (iConfig.getParameter<std::vector<std::string> >("jecPayloadNames") ),
  jecUncName_      (iConfig.getParameter<std::string>("jecUncName"))

{
  jetstoken_          = consumes<std::vector<pat::Jet> >(jets_); 
  globaljetstoken_    = consumes<std::vector<pat::Jet> >(globaljets_);
  genParticlestoken_  = consumes<std::vector<reco::GenParticle> >(genParticles_);
  metToken_           = consumes<std::vector<pat::MET> >(met_);
  globalmetToken_     = consumes<std::vector<pat::MET> >(globalmet_);
  tok_Vertex_         = consumes<reco::VertexCollection>(pvCollection_);  
  tok_beamspot_       = consumes<reco::BeamSpot>(edm::InputTag(bsCollection_));
  pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_));   
  rhoSrc_token_       = consumes<double>(rhoSrc_);

   //Few Counters
  _Nevents_processed = 0.;
  _Nevents_4jets     = 0.;
  _Nevents_4bjets    = 0.;
  _Nevents_ptpass    = 0.;
  _Nevents_mpairs    = 0.;
  _Nevents_deltaM    = 0.;
  _Nevents_passed    = 0.;

  _Noutputs = 0;

  //Few more variables
  is_data        = 1;
  is_not_data    = 0;
  is_mc          = 1;
  is_not_mc      = 0;
  is_json        = 1.0;
  is_silver_json = 1.0;
  is_run         = 1;
  is_xsec        = 1.0;
  is_trigger     = 1.0;
  is_mcFlavour   = 1;
  //nGenJets     = 1;
  jet_id         = 1;
  jet_mcMatchId  =1;


  //Create the histograms and let TFileService handle them
  create_Histos_and_Trees();
}


HAA4bAnalysis::~HAA4bAnalysis()
{
}

// ------------ method called for each event  ------------
void HAA4bAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  _Nevents_processed++;

  bool _show_output(true);
  _Noutputs < 11 ? _Noutputs++ : _show_output = false;

  // define a jet handle and get the jets
  edm::Handle<std::vector<pat::Jet> > jets;  
  iEvent.getByLabel(jets_, jets);
 
  //No point to continue if there aren't 4 jets
   if(jets->size() < 4) return;

  //define a MET handle and get the MET
  edm::Handle<std::vector<pat::MET> > mets;  
  iEvent.getByLabel(met_, mets);

  // Get the mean pt per unit area ("rho")
   edm::Handle<double> rhos;
   iEvent.getByToken(rhoSrc_token_, rhos);

   //Jet Energy Resolution and Scale Factor Stuff
   Resolution    = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
   SF_Resolution = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
   iSetup.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl);

  //Apply Jet Corrections On Fly
  std::vector<JetCorrectorParameters> vPar;
  for ( std::vector<std::string>::const_iterator payloadBegin = jecPayloadNames_.begin(),
        payloadEnd = jecPayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
        JetCorrectorParameters pars(*ipayload);
        vPar.push_back(pars);
    }

  // Make the FactorizedJetCorrector and Uncertainty
   jec_    = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) ); 
   jecUnc_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecUncName_) );
 
   //Alternative way of getting Unceratainity directly from DB without use ot .txt files
   JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
   JetCorrectionUncertainty *jecUnc         = new JetCorrectionUncertainty(JetCorPar);

  //get the Handle of the primary vertex collection and remove the beamspot
  edm::Handle<reco::BeamSpot> bmspot;
  iEvent.getByLabel(bsCollection_,bmspot);

  edm::Handle<reco::VertexCollection> pvcoll;  
  iEvent.getByLabel(pvCollection_,pvcoll);
 
  //require in the event that there is at least one reconstructed vertex
  _nPv =    0;
  npT  = -1.0;
  npIT = -1.0;

  if(pvcoll->size()<=0) return;
  for(reco::VertexCollection::const_iterator vtx=pvcoll->begin();vtx!=pvcoll->end();++vtx) {
    // check that the primary vertex is not a fake one, that is the beamspot (it happens when no primary vertex is reconstructed)
    if(!vtx->isFake()) {
      _nPv++;
    }
  } 


  // PileUp code for examining the Pileup information
  PU_Weight     = 1.;
  PU_WeightUp   = 1;
  PU_WeightDown = 1;

  if (!runningOnData_){
    std::cout<<"Running on Mote Carlo "<<std::endl;
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByLabel(PileupSrc_, PupInfo);  
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	npT  = PVI->getTrueNumInteractions();
	npIT = PVI->getPU_NumInteractions();
      }
    } 


    // calculate weight using above code
    PU_Weight     = Lumiweights_.weight(npT);
    PU_WeightUp   = LumiweightsUp_.weight(npT); 
    PU_WeightDown = LumiweightsDown_.weight(npT);   

    if(_show_output) 
    std::cout<<"PU_Weight for MC is "    <<PU_Weight<<std::endl;
    std::cout<<"PU_WeightUp for MC is "  <<PU_WeightUp<<std::endl;
    std::cout<<"PU_WeightDown for MC is "<<PU_WeightDown<<std::endl;

  }

  if(runningOnData_ && _show_output){
    std::cout<<"Running on Data "<<std::endl;
    std::cout<<"PU_Weight for Data is "<<PU_Weight<<std::endl;
  }

  _Nevents_4jets++;

  int nj1 = -1;
  int nj2 = -1;
  int nj3 = -1;
  int nj4 = -1;
  int nCounter = 0;

  float jetpt1max = -999.;
  float jetpt2max = -999.;
  float jetpt3max = -999.;
  float jetpt4max = -999.;
  
  // loop over jets and determine the ranking
  // based on pt from 1 to 4
  if(_show_output) std::cout<<"Processing event number " << _Nevents_processed << " with equal to or more than 4 Jets" << std::endl;
  for(auto jet = jets->begin(); jet != jets->end(); ++jet){
    float thept = jet->p4().Pt();

    
    //JEC follow Up procedure
    reco::Candidate::LorentzVector uncorrJet;

    // The pat::Jet "knows" if it has been corrected, so here
    //we can "uncorrect" the entire jet to apply the corrections  we want.
    pat::Jet const * pJet = dynamic_cast<pat::Jet const *>( &*jet );
    if ( pJet != 0 ) {
      uncorrJet = pJet->correctedP4(0);
    }
    // Otherwise, if we do not have pat::Jets on input, we just assume
    // the user has not corrected them upstream and use it as raw. 
          else {
                uncorrJet = jet->p4();
                    }

    jec_->setJetEta( uncorrJet.eta() );
    jec_->setJetPt ( uncorrJet.pt() );
    jec_->setJetE  ( uncorrJet.energy() );
    jec_->setJetA  ( jet->jetArea() );
    jec_->setRho   ( *(rhos) );
    //jec_->setRho   ( *(h_rho.product()) );
    //jec_->setNPV   ( h_pv->size() );

    // JEC factor
     corr = jec_->getCorrection();  

    // Now access the uncertainty on the jet energy correction.
    // Pass the corrected jet pt to the "setJetPt" method. 
   
    // Access the "scale up" uncertainty (+1)
    jecUnc_->setJetEta( uncorrJet.eta() );
    jecUnc_->setJetPt( corr * uncorrJet.pt() );
    corrUp = corr * (1 + fabs(jecUnc_->getUncertainty(1)));

   // Access the "scale down" uncertainty (-1)
     jecUnc_->setJetEta( uncorrJet.eta() );
     jecUnc_->setJetPt( corr * uncorrJet.pt() );
     corrDown = corr * ( 1 - fabs(jecUnc_->getUncertainty(-1)) );

    std::cout<<" Jet Energy correction factor = "         <<corr<<std::endl;
    std::cout<<" Jet Enery correction factor ScaleUp = "  <<corrUp<<std::endl;
    std::cout<<" Jet Energy correction factor ScaleDown " <<corrDown<<std::endl;

   // Now plot the value.   
/*  corrPt->Fill( corr * uncorrJet.pt() );
    corrPtUp->Fill( corrUp * uncorrJet.pt() );
    corrPtDown->Fill( corrDown * uncorrJet.pt() );
 */

     //Jet energy Resolution, parameter initialization and value
      Parameters.set(JME::Binning::JetPt, jet->p4().Pt());
      Parameters.set({JME::Binning::JetEta, jet->p4().eta()});
      Parameters.set({JME::Binning::Rho, *rhos});
      //float  smear = getJER(jet->p4().Eta(), 0);

      // Access scale factor and, up and down variation of the scale factor
      Jet_resolution         = Resolution.getResolution(Parameters); 
      Jet_resolutionSF       = SF_Resolution.getScaleFactor(Parameters);
      Jet_resolutionSF_Up    = SF_Resolution.getScaleFactor(Parameters, Variation::UP);
      Jet_resolutionSF_Down  = SF_Resolution.getScaleFactor(Parameters, Variation::DOWN);

     //Printing JER above Info
      std::cout<<"Jet Resolution = "                       <<Jet_resolution<<std::endl;
      std::cout<<"Jet Scale Factor Nominal  = "            <<Jet_resolutionSF<<std::endl;
      std::cout<<"Jet Scale Factor Scale Up Variation = "  <<Jet_resolutionSF_Up<<std::endl;
      std::cout<<"Jet Scale Factor Scale Down variation = "<<Jet_resolutionSF_Down<<std::endl;
   
      //Alternative way to get JEC Uncertainities
      //Retrieving Jet Energy Correction Uncertainities directly from DB without use of text files
      jecUnc->setJetEta(jet->p4().eta());  //Set eta
      jecUnc->setJetPt(jet->p4().Pt());    //Set Pt Here wehave to use the corrected Pt and don't know how. Please checck it 
 
      JEC_Uncertainity = jecUnc->getUncertainty(true);
      ptCor_ScaleUp    = (jet->p4().Pt())*(1+ (1)*JEC_Uncertainity) ;
      ptCor_ScaleDown  = (jet->p4().Pt())*(1+ (-1)*JEC_Uncertainity) ;

      //Printing JEC Uncertainity Info
      std::cout<<"///////////////////////////////////////"<<std::endl;
      std::cout<<"JEC Uncertainity is " <<JEC_Uncertainity<<std::endl;
      std::cout<<"Pt Nominal is "       <<jet->p4().Pt()<<std::endl;
      std::cout<<"Pt Scale Up is "      <<ptCor_ScaleUp<<std::endl;
      std::cout<<"Pt Scale Down is "    <<ptCor_ScaleDown<<std::endl;


    //Jet Btag and other Info
    float thecsv = jet->bDiscriminator(bdiscr_);
    if(_show_output) std::cout<<"Momentum is " << thept << std::endl;

    if(thecsv < minCSV_) continue;

    if(thept > jetpt1max){
      jetpt4max = jetpt3max;
      jetpt3max = jetpt2max;
      jetpt2max = jetpt1max;
      jetpt1max = thept;
      nj4 = nj3;
      nj3 = nj2;
      nj2 = nj1;
      nj1 = nCounter;
    }
    else if(thept > jetpt2max){
      jetpt4max = jetpt3max;
      jetpt3max = jetpt2max;
      jetpt2max = thept;
      nj4 = nj3;
      nj3 = nj2;
      nj2 = nCounter;
    }
    else if(thept > jetpt3max){
      jetpt4max = jetpt3max;
      jetpt3max = thept;
      nj4 = nj3;
      nj3 = nCounter;
    }
    else if(thept > jetpt4max){
      jetpt4max = thept;
      nj4 = nCounter;
    }

    nCounter++;
  }

  if(nj4 < 0) return;
  _Nevents_4bjets++;

  //Fill the tree for QCD background studies

  //Define a jet handle and get the jets without special selections
  edm::Handle<std::vector<pat::Jet> > globaljets;
  iEvent.getByLabel(globaljets_, globaljets);

  //define a handle to access generator level partciles
  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  if(!runningOnData_) iEvent.getByLabel(genParticles_, genParticles);

  //Access MET and store it in the different vectors
  edm::Handle<std::vector<pat::MET> > globalmets;  ///Met
  iEvent.getByLabel(globalmet_, globalmets);

  // Initialization of Background variables
  // Integer initization
  evt    = 0;
  run    = 0;
  lumi   = 0;
  isData = 0;

  HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v            = 0;
  HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v            = 0;
  HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v = 0;
  HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v = 0;
  HLT_HH4bAll = 0;

  //floats initization
 // xsec = 0.;
  //puWeight = 0.;
 // puWeightUp = 0.; 
 // puWeightDown = 0.;
  xsec        = 0.;
  json        = 0;
  json_silver = 0.;
 
  LHE_weights_scale_wgt = 0.; 
  nLHE_weights_pdf = 0.;
  LHE_weights_pdf_id = 0.; 
  LHE_weights_pdf_wgt = 0.;

  //Storing Information in Background Variables
  //run.push_back(iEvent.id().run());  //change run no. also
  run   = is_run;      //Because the run is not considering the actual run but a '1' as input, might be signifying RUN1
  lumi  = iEvent.id().luminosityBlock();
  evt   = iEvent.id().event();

  HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v            = is_trigger; //need proper values
  HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v            = is_trigger;
  HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v = is_trigger;
  HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v = is_trigger; 
  HLT_HH4bAll = is_trigger;

  xsec        = is_xsec;
  json        = is_json;
  json_silver = is_silver_json; 

  if(runningOnData_){ 
    isData       = is_data;          //data is present
    puWeight     = PU_Weight;         
    puWeightUp   = PU_WeightUp;    
    puWeightDown = PU_WeightDown;   
 
   } 
  else{
    isData                = is_not_data;
    puWeight              = PU_Weight;            
    puWeightUp            = PU_WeightUp;
    puWeightDown          = PU_WeightDown;

    LHE_weights_scale_wgt = PU_Weight; //Temprary and need to be changed onwards
    nLHE_weights_pdf      = PU_Weight;
    LHE_weights_pdf_id    = PU_Weight;
    LHE_weights_pdf_wgt   = PU_Weight;

  }
   
  fill_global_Tree(globaljets, genParticles, globalmets);

  //Select the highest csv jets
  pat::Jet jet1 = jets->at(nj1);
  pat::Jet jet2 = jets->at(nj2);
  pat::Jet jet3 = jets->at(nj3);
  pat::Jet jet4 = jets->at(nj4);

 // std::cout<<"jet1FlavourInfobSize = "<<jet3.jetFlavourInfo().getbHadrons().size()<<std::endl;
 // std::cout<<"jet1FlavourInfocSize = "<<jet3.jetFlavourInfo().getcHadrons().size()<<std::endl;

  LorentzVector jet1_4mom = jet1.p4();
  LorentzVector jet2_4mom = jet2.p4();
  LorentzVector jet3_4mom = jet3.p4();
  LorentzVector jet4_4mom = jet4.p4();

   
  if(jet1_4mom.Pt() < minPt_high_ || jet4_4mom.Pt() < minPt_low_) return;
  _Nevents_ptpass++;

  jet1_4mom_tree->SetPxPyPzE(jet1_4mom.Px(),jet1_4mom.Py(),jet1_4mom.Pz(),jet1_4mom.E());
  jet2_4mom_tree->SetPxPyPzE(jet2_4mom.Px(),jet2_4mom.Py(),jet2_4mom.Pz(),jet2_4mom.E());
  jet3_4mom_tree->SetPxPyPzE(jet3_4mom.Px(),jet3_4mom.Py(),jet3_4mom.Pz(),jet3_4mom.E());
  jet4_4mom_tree->SetPxPyPzE(jet4_4mom.Px(),jet4_4mom.Py(),jet4_4mom.Pz(),jet4_4mom.E());

  // Refuse to continue if no combination is above 130 GeV
  if( !check_combinations(jet1_4mom,jet2_4mom,jet3_4mom,jet4_4mom,130.) ) return;
  _Nevents_mpairs++;

  // Get the best pairing possible of two by two bjets, according to invariant mass
  // Convention:
  //1 -> 12 34
  //2 -> 13 24
  //3 -> 14 23
  int combination_flag = get_best_combination(jet1_4mom,jet2_4mom,jet3_4mom,jet4_4mom);

  LorentzVector p_pair1;
  LorentzVector p_pair2;
  if(combination_flag == 1){
    p_pair1 = jet1_4mom + jet2_4mom;
    p_pair2 = jet3_4mom + jet4_4mom;
  }
  else if(combination_flag == 2){
    p_pair1 = jet1_4mom + jet3_4mom;
    p_pair2 = jet2_4mom + jet4_4mom;
  }
  else if(combination_flag == 3){
    p_pair1 = jet1_4mom + jet4_4mom;
    p_pair2 = jet2_4mom + jet3_4mom;
  }

  // avoid a case in which the difference between the masses is above 100 GeV
  if(fabs (p_pair1.M() - p_pair2.M() ) > 100.) return;
  _Nevents_deltaM++;

  // angular distributions between the two pairs
  float var_delta_Phi_pair = p_pair1.Phi() - p_pair2.Phi();
  float var_delta_Eta_pair = p_pair1.Eta() - p_pair2.Eta();

  float total_Mass_4b = (p_pair1 + p_pair2).M();
  if(total_Mass_4b < 260.) return;

  _Nevents_passed++;

  //Make a mass constrained fit
  TKinFitter fitted_Cand = get_fitted_candidate(jet1, jet2, jet3, jet4, combination_flag);

  TLorentzVector jet1_fit = *(fitted_Cand.get4Vec(0));
  TLorentzVector jet2_fit = *(fitted_Cand.get4Vec(1));
  TLorentzVector jet3_fit = *(fitted_Cand.get4Vec(2));
  TLorentzVector jet4_fit = *(fitted_Cand.get4Vec(3));
  
  float_t total_Mass_4b_fitted = (jet1_fit + jet2_fit + jet3_fit + jet4_fit).M();

  std::cout << "4b mass before fit = " << total_Mass_4b << std::endl;
  std::cout << "4b mass after  fit = " << total_Mass_4b_fitted << std::endl;

  jet1_4mom_tree_fit->SetPxPyPzE(jet1_fit.Px(),jet1_fit.Py(),jet1_fit.Pz(),jet1_fit.E());
  jet2_4mom_tree_fit->SetPxPyPzE(jet2_fit.Px(),jet2_fit.Py(),jet2_fit.Pz(),jet2_fit.E());
  jet3_4mom_tree_fit->SetPxPyPzE(jet3_fit.Px(),jet3_fit.Py(),jet3_fit.Pz(),jet3_fit.E());
  jet4_4mom_tree_fit->SetPxPyPzE(jet4_fit.Px(),jet4_fit.Py(),jet4_fit.Pz(),jet4_fit.E());

  // now plot a few quantities for the jets
  var_jet1Btag = jet1.bDiscriminator(bdiscr_);
  var_jet2Btag = jet2.bDiscriminator(bdiscr_);
  var_jet3Btag = jet3.bDiscriminator(bdiscr_);
  var_jet4Btag = jet4.bDiscriminator(bdiscr_);

  h_jet1pt->Fill(jet1.p4().Pt());
  h_jet2pt->Fill(jet2.p4().Pt());
  h_jet3pt->Fill(jet3.p4().Pt());
  h_jet4pt->Fill(jet4.p4().Pt());

  h_jet1eta->Fill(jet1.p4().Eta());
  h_jet2eta->Fill(jet2.p4().Eta());
  h_jet3eta->Fill(jet3.p4().Eta());
  h_jet4eta->Fill(jet4.p4().Eta());

  h_jet1Btag->Fill(var_jet1Btag);
  h_jet2Btag->Fill(var_jet2Btag);
  h_jet3Btag->Fill(var_jet3Btag);
  h_jet4Btag->Fill(var_jet4Btag);

  // and the invariant masses
  h_m_pair1->Fill(p_pair1.M());
  h_m_pair2->Fill(p_pair2.M());
  h_m_4b->Fill(total_Mass_4b);
  h_m_4b_fitted->Fill(total_Mass_4b_fitted);
  h_m4b_m12->Fill(total_Mass_4b,p_pair1.M());
  h_m4b_m34->Fill(total_Mass_4b,p_pair2.M());

  // and the angular "distributions"
  h_delta_Phi_pair->Fill(var_delta_Phi_pair);
  h_delta_Eta_pair->Fill(var_delta_Eta_pair);
 
  // and the PU info
  h_nPv->Fill(_nPv); 
  h_Rw_nPv->Fill(float(_nPv)-1,float(PU_Weight)); // subtract primary interaction
  h_PUTrue->Fill(npT);
  h_PUInTime->Fill(npIT);
  
  // once you have the event weight, you can plot reweighted distributions of important event quantities 
  h_PUWeight->Fill(PU_Weight);
  h_Rw_PUTrue->Fill(npT, PU_Weight);
  h_Rw_PUInTime->Fill(npIT, PU_Weight);
  
  mytree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void HAA4bAnalysis::beginJob()
{
  // Flag for PileUp reweighting
  if (!runningOnData_)
   Lumiweights_     = edm::LumiReWeighting("pileUpHistogramFromjson_Nominal.root",   "MCpileUp_25ns_Recent2016.root", "pileup", "pileup");
   LumiweightsUp_   = edm::LumiReWeighting("pileUpHistogramFromjson_ScaleUp.root",   "MCpileUp_25ns_Recent2016.root", "pileup", "pileup");
   LumiweightsDown_ = edm::LumiReWeighting("pileUpHistogramFromjson_ScaleDown.root", "MCpileUp_25ns_Recent2016.root", "pileup", "pileup");
 
}

// ------------ method called once each job just after ending the event loop  ------------
void HAA4bAnalysis::endJob() 
{
  h_Events->Fill(0.5,_Nevents_processed);
  h_Events->Fill(1.5,_Nevents_4jets);
  h_Events->Fill(2.5,_Nevents_4bjets);
  h_Events->Fill(3.5,_Nevents_ptpass);
  h_Events->Fill(4.5,_Nevents_mpairs);
  h_Events->Fill(5.5,_Nevents_deltaM);
  h_Events->Fill(6.5,_Nevents_passed);

  std::cout << "###########################################"      << std::endl;
  std::cout << "EVENT SUMMARY"                                    << std::endl;
  std::cout << "Number of events processed = "                    << _Nevents_processed << std::endl;
  std::cout << "Number of events with 4 jets = "                  << _Nevents_4jets << std::endl;
  std::cout << "Number of events with 3 bjets = "                 << _Nevents_4bjets << std::endl;
  std::cout << "Number of events passing pt cuts = "              << _Nevents_ptpass << std::endl;
  std::cout << "Number of events with m_12/34 > 120 GeV = "       << _Nevents_mpairs << std::endl;
  std::cout << "Number of events with delta_m_12/34 < 100 GeV = " << _Nevents_deltaM << std::endl;
  std::cout << "Number of events passing the selection = "        << _Nevents_passed << std::endl;
  std::cout << "#######################"                          << std::endl;
}

int HAA4bAnalysis::get_best_combination(LorentzVector& m1, LorentzVector& m2, LorentzVector& m3, LorentzVector& m4){

  float diff_m_12_34 = fabs( (m1+m2).M() - (m3+m4).M() );
  float diff_m_13_24 = fabs( (m1+m3).M() - (m2+m4).M() );
  float diff_m_14_23 = fabs( (m1+m4).M() - (m2+m3).M() );

  std::vector<float> diff_vector;

  diff_vector.push_back(diff_m_12_34);
  diff_vector.push_back(diff_m_13_24);
  diff_vector.push_back(diff_m_14_23);
  std::sort(diff_vector.begin(), diff_vector.end());

  if(diff_vector.at(0) == diff_m_12_34) return 1;
  else if(diff_vector.at(0) == diff_m_13_24) return 2;
  else if(diff_vector.at(0) == diff_m_14_23) return 3;
  else assert(0);

  return -1;
}

bool HAA4bAnalysis::check_combinations(LorentzVector& m1, LorentzVector& m2, LorentzVector& m3, LorentzVector& m4, float mcut){

  bool b12 = (m1+m2).M() > mcut;
  bool b13 = (m1+m3).M() > mcut;
  bool b14 = (m1+m4).M() > mcut;
  bool b34 = (m3+m4).M() > mcut;
  bool b24 = (m2+m4).M() > mcut;
  bool b23 = (m2+m3).M() > mcut;

  return b12 || b13 || b14 || b34 || b24 || b23;
}

TKinFitter HAA4bAnalysis::get_fitted_candidate(pat::Jet& Jet1, pat::Jet& Jet2, pat::Jet& Jet3, pat::Jet& Jet4, int best_combination)
{
  static int debug_counter = 0;

  TLorentzVector jet1_4mom(Jet1.px(),Jet1.py(),Jet1.pz(),Jet1.energy()) ;
  TLorentzVector jet2_4mom(Jet2.px(),Jet2.py(),Jet2.pz(),Jet2.energy()) ;
  TLorentzVector jet3_4mom(Jet3.px(),Jet3.py(),Jet3.pz(),Jet3.energy()) ;
  TLorentzVector jet4_4mom(Jet4.px(),Jet4.py(),Jet4.pz(),Jet4.energy()) ;

  //Covariance matrix on get uncertainties
  TMatrixD m1(3,3);
  TMatrixD m2(3,3);
  TMatrixD m3(3,3);
  TMatrixD m4(3,3);
  m1.Zero();
  m2.Zero();
  m3.Zero();
  m4.Zero();

  //FIXME: for now, just put 5% to test, but we need to put the proper jet energy resolution
  m1(0,0) = 0.08 * Jet1.p4().Et();
  m1(1,1) = 0.01 * Jet1.p4().Eta();
  m1(2,2) = 0.01 * Jet1.p4().Phi();

  m2(0,0) = 0.08 * Jet2.p4().Et();
  m2(1,1) = 0.01 * Jet2.p4().Eta();
  m2(2,2) = 0.01 * Jet2.p4().Phi();

  m3(0,0) = 0.08 * Jet3.p4().Et();
  m3(1,1) = 0.01 * Jet3.p4().Eta();
  m3(2,2) = 0.01 * Jet3.p4().Phi();

  m4(0,0) = 0.08 * Jet4.p4().Et();
  m4(1,1) = 0.01 * Jet4.p4().Eta();
  m4(2,2) = 0.01 * Jet4.p4().Phi();

  if(debug_counter < 5){
    std::cout << "Jet1 resolution pt  = " << m1(0,0) << std::endl;
    std::cout << "Jet1 resolution eta = " << m1(1,1) << std::endl;
    std::cout << "Jet1 resolution phi = " << m1(2,2) << std::endl;
  }

  TFitParticleEtEtaPhi jet1( "Jet1", "Jet1", &jet1_4mom, &m1 );
  TFitParticleEtEtaPhi jet2( "Jet2", "Jet2", &jet2_4mom, &m2 );
  TFitParticleEtEtaPhi jet3( "Jet3", "Jet3", &jet3_4mom, &m3 );
  TFitParticleEtEtaPhi jet4( "Jet4", "Jet4", &jet4_4mom, &m4 );

  TLorentzVector p_pair1;
  TLorentzVector p_pair2;
  if(best_combination == 1){
    p_pair1 = jet1_4mom + jet2_4mom;
    p_pair2 = jet3_4mom + jet4_4mom;
  }
  else if(best_combination == 2){
    p_pair1 = jet1_4mom + jet3_4mom;
    p_pair2 = jet2_4mom + jet4_4mom;
  }
  else if(best_combination == 3){
    p_pair1 = jet1_4mom + jet4_4mom;
    p_pair2 = jet2_4mom + jet3_4mom;
  }

  //Best estimate of A->bb is given by the average of the two
  float mA_estimate = (p_pair1.M()+p_pair2.M())/2.;
  if(debug_counter < 5) std::cout << "mA estimate = " << mA_estimate << std::endl;

  //Create the constraints
  TFitConstraintM mCons1( "A1MassConstraint", "A1Mass-Constraint", 0, 0 , 0.); //request to have mass difference = 0

  if(best_combination == 1){
    mCons1.addParticles1( &jet1, &jet2 );
    mCons1.addParticles2( &jet3, &jet4 );
  }
  else if(best_combination == 2){
    mCons1.addParticles1( &jet1, &jet3 );
    mCons1.addParticles2( &jet2, &jet4 );
  }
  else if(best_combination == 3){
    mCons1.addParticles1( &jet1, &jet4 );
    mCons1.addParticles2( &jet2, &jet3 );
  }

  //Definition of the fitter
  //Add four measured particles(jets)
  //Add two constraints
  TKinFitter fitter("fitter", "fitter");
  fitter.addMeasParticle( &jet1  );
  fitter.addMeasParticle( &jet2  );
  fitter.addMeasParticle( &jet3  );
  fitter.addMeasParticle( &jet4  );
  fitter.addConstraint  ( &mCons1);

  //Set convergence criteria
  fitter.setMaxNbIter( 40 );
  fitter.setMaxDeltaS( 1e-2 );
  fitter.setMaxF( 1e-1 );
  if(debug_counter < 5) fitter.setVerbosity(1);

  //Perform the fit
  if(debug_counter < 5) std::cout << "Performing kinematic fit..." << std::endl;
  fitter.fit();
  if(debug_counter < 5) std::cout << "Done." << std::endl;

  if(debug_counter < 5){
    std::cout << "===================================" << std ::endl;
    std::cout << "-> Number of measured Particles  : " << fitter.nbMeasParticles() << std::endl;
    std::cout << "-> Number of unmeasured particles: " << fitter.nbUnmeasParticles() << std::endl;
    std::cout << "-> Number of constraints         : " << fitter.nbConstraints() << std::endl;
    std::cout << "-> Number of degrees of freedom  : " << fitter.getNDF() << std::endl;
    std::cout << "-> Number of parameters A        : " << fitter.getNParA() << std::endl;
    std::cout << "-> Number of parameters B        : " << fitter.getNParB() << std::endl;
    std::cout << "-> Maximum number of iterations  : " << fitter.getMaxNumberIter() << std::endl;
    std::cout << "-> Maximum deltaS                : " << fitter.getMaxDeltaS() << std::endl;
    std::cout << "-> Maximum F                     : " << fitter.getMaxF() << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++" << std ::endl;
    std::cout << "-> Status                        : " << fitter.getStatus() << std::endl;
    std::cout << "-> Number of iterations          : " << fitter.getNbIter() << std::endl;
    std::cout << "-> S                             : " << fitter.getS() << std::endl;
    std::cout << "-> F                             : " << fitter.getF() << std::endl;
    std::cout << "===================================" << std ::endl;
  }

  debug_counter++;

  return fitter;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HAA4bAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void HAA4bAnalysis::fill_global_Tree(edm::Handle<std::vector<pat::Jet> >& globaljets, edm::Handle<std::vector<reco::GenParticle> >& genParticles, edm::Handle<std::vector<pat::MET> > &globalmets){


  // setup calibration + reader //setup calibration readers (once)   
   BTagCalibration calib("csvv2", "CSVv2_Moriond17_B_H.csv");
   BTagCalibration calib_cs("csvv2", "CSVv2_Moriond17_B_H.csv");

   BTagCalibrationReader reader(BTagEntry::OP_MEDIUM,        // operating point
                               "central",                   // central sys type
                               {"up", "down"});             // other sys types

   reader.load(calib,                                       // calibration instance
               BTagEntry::FLAV_B,                           // btag flavour
              "comb");                                      // measurement type
  
   //Jet energy scale, "jes"; JESUp and JESDown
   BTagCalibrationReader reader_JES(BTagEntry::OP_RESHAPING,"central",{"up_jes", "down_jes"});  
   reader_JES.load(calib_cs, BTagEntry::FLAV_B, "iterativefit");    

   //Light flavor contamination, "lf"; Light Flavour and Scale variation
   BTagCalibrationReader reader_LF(BTagEntry::OP_RESHAPING,"central",{"up_lf", "down_lf"});
   reader_LF.load(calib_cs, BTagEntry::FLAV_B, "iterativefit");

   //Heavy flavor contamination, "hf"; Heavy Flavour and Scale variation
   BTagCalibrationReader reader_HF(BTagEntry::OP_RESHAPING,"central",{"up_hf", "down_hf"});  
   reader_HF.load(calib_cs, BTagEntry::FLAV_UDSG, "iterativefit");

   //Linear statistical fluctuations, "hfstats1"; HFStats1Up and  HFStats1Down
   BTagCalibrationReader reader_HFStats1(BTagEntry::OP_RESHAPING,"central",{"up_hfstats1", "down_hfstats1"});    
   reader_HFStats1.load(calib_cs, BTagEntry::FLAV_B, "iterativefit");

   //Quadratic statistical fluctuations, "hfstats2"; HFStats2Up and  HFStats2Down
   BTagCalibrationReader reader_HFStats2(BTagEntry::OP_RESHAPING,"central",{"up_hfstats2", "down_hfstats2"}); 
   reader_HFStats2.load(calib_cs, BTagEntry::FLAV_B, "iterativefit");

   //Linear statistical fluctuations, "lfstats1"; LFStats1Up and  LFStats1Down
   BTagCalibrationReader reader_LFStats1(BTagEntry::OP_RESHAPING,"central",{"up_lfstats1", "down_lfstats1"});
   reader_LFStats1.load(calib_cs, BTagEntry::FLAV_UDSG, "iterativefit");


   //Quadratic statistical fluctuations, "lfstats2"; HFStats2Up and  HFStats2Down
   BTagCalibrationReader reader_LFStats2(BTagEntry::OP_RESHAPING,"central",{"up_lfstats2", "down_lfstats2"});
   reader_LFStats2.load(calib_cs, BTagEntry::FLAV_UDSG, "iterativefit");
 
   //Linear uncertainties, "cferr1" for charm flavour; CFErr1Up and CFErr1Down
   BTagCalibrationReader reader_CFErr1(BTagEntry::OP_RESHAPING,"central",{"up_cferr1", "down_cferr1"});
   reader_CFErr1.load(calib_cs, BTagEntry::FLAV_C, "iterativefit");

   //Quadratic uncertainties, "cferr2"; CFErr2Up and CFErr2Down
   BTagCalibrationReader reader_CFErr2(BTagEntry::OP_RESHAPING,"central",{"up_cferr2", "down_cferr2"}); 
   reader_CFErr2.load(calib_cs, BTagEntry::FLAV_C, "iterativefit");


  //Few more background variables and their initialization
  Jet_pt   = 0.;
  Jet_eta  = 0.;
  Jet_phi  = 0.;
  Jet_mass = 0.;
  Jet_btag = 0.;  

  bTagWeight_LFUp = 0.;
  bTagWeight_LFDown = 0.;
  bTagWeight_LFStats1Up = 0.;
  bTagWeight_LFStats1Down = 0.;
  bTagWeight_LFStats2Down = 0.;
  bTagWeight_LFStats2Up = 0.;

  bTagWeight_HFUp = 0.;
  bTagWeight_HFDown = 0.;
  bTagWeight_HFStats1Up = 0.;
  bTagWeight_HFStats1Down = 0.;
  bTagWeight_HFStats2Up = 0.;
  bTagWeight_HFStats2Down = 0.;

  bTagWeight_cErr1Down = 0.;
  bTagWeight_cErr1Up = 0.;
  bTagWeight_cErr2Up = 0.;
  bTagWeight_cErr2Down = 0.;

  bTagWeight_JESDown = 0.;
  bTagWeight_JESUp   = 0.;
  bTagWeight         = 0.;

  Jet_bTagWeight = 0.;

  Jet_btagCSV = 0. ;
  Jet_btagCMVA = 0.; 
  Jet_btagCSVV0 = 0.; 
  Jet_btagCMVAV2 = 0.; 
  Jet_mcPt = 0. ;
 
  Jet_corr         = 0.; 
  Jet_corr_JECUp   = 0.; 
  Jet_corr_JECDown = 0.; 

  Jet_corr_JER     = 0.;
  Jet_corr_JERUp   = 0.; 
  Jet_corr_JERDown = 0.;


  for(auto jet = globaljets->begin(); jet != globaljets->end(); ++jet){
    float thept   = jet->p4().Pt();
    float theeta  = jet->p4().Eta();
    float thephi  = jet->p4().Phi();
    float themass = jet->p4().M();
    float thecsv  = jet->bDiscriminator(bdiscr_);

    //Jet Hadron and Parton Flavour
    Jet_partonFlavour = jet->hadronFlavour();
    Jet_hadronFlavour = jet->partonFlavour();

    //Store info into variables and also need to control the loop as the data might be over written during the next pass of the above loop
    Jet_pt   = thept;
    Jet_eta  = theeta;
    Jet_phi  = thephi;
    Jet_mass = themass;
    Jet_btag = thecsv;

   if(!runningOnData_){

     //Btag Scale factors and scale Variations
     //Eta values outside -2.4 < eta < 2.4 and discriminator values outside 0 < discr < 1, Pt <20 will yield a scale factor of 0
     //using eval_auto_bounds() function
     bjet_scaleFactor           = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().Pt()); 
     bjet_scaleFactor_ScaleUp   = reader.eval_auto_bounds("up",      BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().Pt());
     bjet_scaleFactor_ScaleDown = reader.eval_auto_bounds("down",    BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt()); 

     std::cout<<"Btaging Scale factors"       <<std::endl;
     std::cout<<"Btag Scale factor central "  <<bjet_scaleFactor <<std::endl;
     std::cout<<"Btag Scale factor scaleUp "  <<bjet_scaleFactor_ScaleUp<<std::endl;
     std::cout<<"Btag Scale factor scaleDown" <<bjet_scaleFactor_ScaleDown<<std::endl;

    //Jet Energy Scale Factors and Scale variation
    bTagWeight_JESUp   = reader_JES.eval_auto_bounds("up_jes",   BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());
    bTagWeight_JESDown = reader_JES.eval_auto_bounds("down_jes", BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"Jet Energy ScaleUp factor "    <<bTagWeight_JESUp<<std::endl;
    std::cout<<"Jet Energy ScaleDown factor "  <<bTagWeight_JESDown<<std::endl;

    //Light Flavour and Scale Variations
    Jet_bTagWeightLFUp   = reader_LF.eval_auto_bounds("up_lf",   BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightLFDown = reader_LF.eval_auto_bounds("down_lf", BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());
 
    std::cout<<"Light Flavour ScaleUp factor  "  <<Jet_bTagWeightLFUp<<std::endl;
    std::cout<<"Light Flavour ScaleDown factor " <<Jet_bTagWeightLFDown<<std::endl;

    //Heavy Flavour and Scale Variations
    Jet_bTagWeightHFUp   = reader_HF.eval_auto_bounds("up_hf",   BTagEntry::FLAV_UDSG, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightHFDown = reader_HF.eval_auto_bounds("down_hf", BTagEntry::FLAV_UDSG, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"Heavy Flavour ScaleUp factor  "      <<Jet_bTagWeightHFUp<<std::endl;
    std::cout<<"Heavy Flavour ScaleDown factor "     <<Jet_bTagWeightHFDown<<std::endl;

    //// HFStats1; Up and down variations
    Jet_bTagWeightHFStats1Up   = reader_HFStats1.eval_auto_bounds("up_hfstats1",   BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightHFStats1Down = reader_HFStats1.eval_auto_bounds("down_hfstats1", BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"Heavy Flavour hfstats1 ScaleUp factor  "      <<Jet_bTagWeightHFStats1Up<<std::endl;
    std::cout<<"Heavy Flavour hfstats1 ScaleDown factor "     <<Jet_bTagWeightHFStats1Down<<std::endl;

    //// HFStats2; Up and down variations
    Jet_bTagWeightHFStats2Up   = reader_HFStats2.eval_auto_bounds("up_hfstats2",   BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightHFStats2Down = reader_HFStats2.eval_auto_bounds("down_hfstats2", BTagEntry::FLAV_B, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"Heavy Flavour hfstats2 ScaleUp factor  "      <<Jet_bTagWeightHFStats2Up<<std::endl;
    std::cout<<"Heavy Flavour hfstats2 ScaleDown factor "     <<Jet_bTagWeightHFStats2Down<<std::endl;

    // LFStats1; Up and down variations
    Jet_bTagWeightLFStats1Up   = reader_LFStats1.eval_auto_bounds("up_lfstats1",   BTagEntry::FLAV_UDSG, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightLFStats1Down = reader_LFStats1.eval_auto_bounds("down_lfstats1", BTagEntry::FLAV_UDSG, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"Light Flavour lfstats1 ScaleUp factor  "      <<Jet_bTagWeightLFStats1Up<<std::endl;
    std::cout<<"Light Flavour lfstats1 ScaleDown factor "     <<Jet_bTagWeightLFStats1Down<<std::endl;

    // LFStats2; Up and down variations
    Jet_bTagWeightLFStats2Up   = reader_LFStats2.eval_auto_bounds("up_lfstats2",   BTagEntry::FLAV_UDSG, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightLFStats2Down = reader_LFStats2.eval_auto_bounds("down_lfstats2", BTagEntry::FLAV_UDSG, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"Light Flavour lfstats2 ScaleUp factor  "      <<Jet_bTagWeightLFStats2Up<<std::endl;
    std::cout<<"Light Flavour lfstats2 ScaleDown factor "     <<Jet_bTagWeightLFStats2Down<<std::endl;

    //CFErr1Up and CFErr1Down
    Jet_bTagWeightcErr1Up   = reader_CFErr1.eval_auto_bounds("up_cferr1",   BTagEntry::FLAV_C, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightcErr1Down = reader_CFErr1.eval_auto_bounds("down_cferr1", BTagEntry::FLAV_C, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"cferr1 ScaleUp factor  "      <<Jet_bTagWeightcErr1Up<<std::endl;
    std::cout<<"cferr1 ScaleDown factor "     <<Jet_bTagWeightcErr1Down<<std::endl;

    //CFErr2Up and CFErr2Down
    Jet_bTagWeightcErr2Up   = reader_CFErr2.eval_auto_bounds("up_cferr2",   BTagEntry::FLAV_C, jet->p4().Eta(), jet->p4().pt());
    Jet_bTagWeightcErr2Down = reader_CFErr2.eval_auto_bounds("down_cferr2", BTagEntry::FLAV_C, jet->p4().Eta(), jet->p4().pt());

    std::cout<<"cferr2 ScaleUp factor  "      <<Jet_bTagWeightcErr2Up <<std::endl;
    std::cout<<"cferr2 ScaleDown factor "     <<Jet_bTagWeightcErr2Down<<std::endl;

    bTagWeight_LFUp = thecsv;  //Need to put a proper value onwards and also to differentiate between bTagWeight and Jet_bTagWeight etc
    bTagWeight_LFDown = thecsv;
    bTagWeight_LFStats1Up = thecsv;
    bTagWeight_LFStats1Down = thecsv;
    bTagWeight_LFStats2Down = thecsv;
    bTagWeight_LFStats2Up = thecsv;

    bTagWeight_HFUp = thecsv;
    bTagWeight_HFDown = thecsv;
    bTagWeight_HFStats1Up = thecsv;
    bTagWeight_HFStats1Down = thecsv;
    bTagWeight_HFStats2Up = thecsv;
    bTagWeight_HFStats2Down = thecsv;

    bTagWeight_cErr1Down = thecsv;
    bTagWeight_cErr1Up = thecsv;
    bTagWeight_cErr2Up = thecsv;
    bTagWeight_cErr2Down = thecsv;

    bTagWeight_JESDown = Jet_resolutionSF_Down; //has to be changed
    bTagWeight_JESUp   = Jet_resolutionSF_Up ;

    bTagWeight = bjet_scaleFactor;

      } // Extend this braket to include more variables given below if they are also to be applied on MC


    Jet_btagCSV = thecsv;
    Jet_btagCMVA = thecsv; 
    Jet_btagCSVV0 = thecsv; 
    Jet_btagCMVAV2 = thecsv; 
    Jet_mcPt = thecsv;

    Jet_corr          = corr; 
    Jet_corr_JECUp    = corrUp; 
    Jet_corr_JECDown  = corrDown; 

    Jet_corr_JER      = Jet_resolutionSF;     //Please check whether we have to use Scale factor or resolution
    Jet_corr_JERUp    = Jet_resolutionSF_Up;  //Similarly here
    Jet_corr_JERDown  = Jet_resolutionSF_Down;

  }
  //Generator level info from MC
  genWeight   = 0.;

  GenJet_pt   = 0.; 
  GenJet_eta  = 0.; 
  GenJet_phi  = 0.; 
  GenJet_mass = 0.;

  Jet_id = 0;
  Jet_mcFlavour = 0;
  nGenJet = 0; //number of Gen Jets
  Jet_mcMatchId = 0;

  //Gen b-infos
  Genb_pt.clear();
  Genb_eta.clear();
  Genb_phi.clear();
  Genb_mass.clear();

  if(!runningOnData_){
    nGenJets = 0;
    for(auto genpart = genParticles->begin(); genpart != genParticles->end(); genpart++){
       
      nGenJets++;
      genWeight     = is_mc; // will have to fix the loop problem as the data will get over written for next pass
      GenJet_pt     = genpart->pt(); 
      GenJet_eta    = genpart->eta(); 
      GenJet_phi    = genpart->phi(); 
      GenJet_mass   = genpart->mass();
      Jet_id        = jet_id;
      Jet_mcFlavour = is_mcFlavour;  
      nGenJet       = nGenJets;  //Need to look for it if it is correct
      Jet_mcMatchId = jet_mcMatchId;   // same as Ref 3

      // Signal Information from MC
      if( abs(genpart->pdgId()) != 5) continue;
      if( abs(genpart->mother(0)->pdgId()) != 36) continue; //A -> bb

      Genb_pt.push_back(genpart->pt());
      Genb_eta.push_back(genpart->eta());
      Genb_phi.push_back(genpart->phi());
      Genb_mass.push_back(genpart->mass());

    }

    // MET information for background
    met_pt   = 0.;
    met_eta  = 0.; 
    met_phi  = 0.; 
    met_mass = 0.;
 
    float mthept   = (globalmets->front()).pt();  //mthept-->Met the pt (pt of met)
    float mtheeta  = (globalmets->front()).eta();
    float mthephi  = (globalmets->front()).phi();
    float mthemass = (globalmets->front()).mass();
   
    met_pt   = mthept;
    met_eta  = mtheeta; 
    met_phi  = mthephi; 
    met_mass = mthemass;
   
  }

  //Fill global tree
  globalTree->Fill();

}

void HAA4bAnalysis::create_Histos_and_Trees(){
  h_jet1pt = fs->make<TH1F>("h_jet1pt", "P_t of the 1st jet", 200, 0., 500.);
  h_jet2pt = fs->make<TH1F>("h_jet2pt", "P_t of the 2nd jet", 200, 0., 500.);
  h_jet3pt = fs->make<TH1F>("h_jet3pt", "P_t of the 3rd jet", 200, 0., 500.);
  h_jet4pt = fs->make<TH1F>("h_jet4pt", "P_t of the 4th jet", 200, 0., 500.);

  h_jet1eta = fs->make<TH1F>("h_jet1eta", "#eta of the 1st jet", 50, -10., 10.);
  h_jet2eta = fs->make<TH1F>("h_jet2eta", "#eta of the 2nd jet", 50, -10., 10.);
  h_jet3eta = fs->make<TH1F>("h_jet3eta", "#eta of the 3rd jet", 50, -10., 10.);
  h_jet4eta = fs->make<TH1F>("h_jet4eta", "#eta of the 4th jet", 50, -10., 10.);

  h_jet1Btag = fs->make<TH1F>("h_jet1Btag", "B-tag discriminant of the 1st jet", 50, minCSV_, 1.);
  h_jet2Btag = fs->make<TH1F>("h_jet2Btag", "B-tag discriminant of the 2nd jet", 50, minCSV_, 1.);
  h_jet3Btag = fs->make<TH1F>("h_jet3Btag", "B-tag discriminant of the 3rd jet", 50, minCSV_, 1.);
  h_jet4Btag = fs->make<TH1F>("h_jet4Btag", "B-tag discriminant of the 4th jet", 50, minCSV_, 1.);

  h_m_pair1     = fs->make<TH1F>("h_m_pair1", "Invariant mass of the first jet pair", 200, 130., 800.);
  h_m_pair2     = fs->make<TH1F>("h_m_pair2", "Invariant mass of the second jet pair", 200, 130., 800.);
  h_m_4b        = fs->make<TH1F>("h_m_4b", "Invariant mass of the four b jets", 200, 130., 1000.);
  h_m_4b_fitted = fs->make<TH1F>("h_m_4b_fitted", "Invariant mass of the four b jets after fit", 200, 130., 1000.);
  h_m4b_m12     = fs->make<TH2F>("h_m4b_m12", "Invariant mass of 4b vs m12", 200, 130., 1000.,200,130.,800.);
  h_m4b_m34     = fs->make<TH2F>("h_m4b_m34", "Invariant mass of 4b vs m34", 200, 130., 1000.,200,130.,800.);

  h_delta_Phi_pair = fs->make<TH1F>("h_delta_Phi_pair", "#Delta_{#phi} between the two jet pairs", 50, 0., 6.28);
  h_delta_Eta_pair = fs->make<TH1F>("h_delta_Eta_pair", "#Delta_{#eta} between the two jet pairs", 50, -10., 10.);

  h_Events = fs->make<TH1F>("h_Events", "Event counting in different steps", 7, 0., 7.);

  h_nPv         = fs->make<TH1F>("h_nPv", "No. of primary verticies", 50, 0., 50.);              //Few more Histogram inclusions
  h_Rw_nPv      = fs->make<TH1F>("h_Rw_nPv", "Reweighted No. of primary verticies", 50, 0., 50.);
  h_PUInTime    = fs->make<TH1F>("h_PUInTime","Input No. in-time pileup interactions",50,0.,50.);
  h_PUTrue      = fs->make<TH1F>("h_PUTrue","Input True pileup interactions",50,0.,50.);
  h_Rw_PUTrue   = fs->make<TH1F>("h_Rw_PUTrue","Reweighted True pileup interactions",50,0.,50.);
  h_Rw_PUInTime = fs->make<TH1F>("h_Rw_PUInTime","Reweighted in-time pileup interactions",50,0.,50.);
  h_PUWeight    = fs->make<TH1F>("h_PUWeight","Event weight",50,0.,50.);

  jet1_4mom_tree = new TLorentzVector();
  jet2_4mom_tree = new TLorentzVector();
  jet3_4mom_tree = new TLorentzVector();
  jet4_4mom_tree = new TLorentzVector();

  jet1_4mom_tree_fit = new TLorentzVector();
  jet2_4mom_tree_fit = new TLorentzVector();
  jet3_4mom_tree_fit = new TLorentzVector();
  jet4_4mom_tree_fit = new TLorentzVector();

  // create the tree and let TFileService handle it
  mytree = fs->make<TTree>("mytree", "Tree containing events after presel");
  mytree->Branch("jet1_4mom","TLorentzVector",&jet1_4mom_tree);
  mytree->Branch("jet2_4mom","TLorentzVector",&jet2_4mom_tree);
  mytree->Branch("jet3_4mom","TLorentzVector",&jet3_4mom_tree);
  mytree->Branch("jet4_4mom","TLorentzVector",&jet4_4mom_tree);

  mytree->Branch("jet1_4mom_fit","TLorentzVector",&jet1_4mom_tree_fit);
  mytree->Branch("jet2_4mom_fit","TLorentzVector",&jet2_4mom_tree_fit);
  mytree->Branch("jet3_4mom_fit","TLorentzVector",&jet3_4mom_tree_fit);
  mytree->Branch("jet4_4mom_fit","TLorentzVector",&jet4_4mom_tree_fit);

  mytree->Branch("jet1Btag",&var_jet1Btag,"jet1Btag/F");
  mytree->Branch("jet2Btag",&var_jet2Btag,"jet2Btag/F");
  mytree->Branch("jet3Btag",&var_jet3Btag,"jet3Btag/F");
  mytree->Branch("jet4Btag",&var_jet4Btag,"jet4Btag/F");

  mytree->Branch("N_nPv", &_nPv, "_nPv/I");                                          //Filling Primary Verticies 


  mytree->Branch("Jet_corr",         &Jet_corr,         "Jet_corr/F");
  mytree->Branch("Jet_corr_JECUp",   &Jet_corr_JECUp,   "Jet_corr_JECUp/F");
  mytree->Branch("Jet_corr_JECDown", &Jet_corr_JECDown, "Jet_corr_JECDown/F");

  mytree->Branch("Jet_corr_JER",     &Jet_corr_JER,     "Jet_corr_JER");
  mytree->Branch("Jet_corr_JERUp",   &Jet_corr_JERUp,   "Jet_corr_JERUp/F");
  mytree->Branch("Jet_corr_JERDown", &Jet_corr_JERDown, "Jet_corr_JERDown/F");
 

  // pileUp histograms
  mytree->Branch("PUTrue",   &npT, "npT/F");
  mytree->Branch("PUInTime", &npIT, "npIT/F");
  mytree->Branch("PUWeight", &PU_Weight, "PU_Weight/F");
  
  globalTree = fs->make<TTree>("globalTree", "Tree containing most of the event info to study background");

  //jet info
  globalTree->Branch("run",  &run,  "run/i" );
  globalTree->Branch("evt",  &evt,  "evt/l");
  globalTree->Branch("lumi", &lumi, "lumi/i");

  globalTree->Branch("Jet_pt",   &Jet_pt,   "Jet_pt/F");
  globalTree->Branch("Jet_eta",  &Jet_eta,  "Jet_eta/F");
  globalTree->Branch("Jet_phi",  &Jet_phi,  "Jet_phi/F");
  globalTree->Branch("Jet_mass", &Jet_mass, "Jet_mass/F");

  globalTree->Branch("Jet_bTagWeight",        &Jet_bTagWeight,        "Jet_bTagWeight/F");
  globalTree->Branch("Jet_bTagWeightJESUp",   &Jet_bTagWeightJESUp,   "Jet_bTagWeightJESUp/F");
  globalTree->Branch("Jet_bTagWeightJESDown", &Jet_bTagWeightJESDown, "Jet_bTagWeightJESDown/F");

  globalTree->Branch("Jet_bTagWeightLFUp",         &Jet_bTagWeightLFUp,         "Jet_bTagWeightLFUp/F");
  globalTree->Branch("Jet_bTagWeightLFDown",       &Jet_bTagWeightLFDown,       "Jet_bTagWeightLFDown/F");
  globalTree->Branch("Jet_bTagWeightLFStats1Up",   &Jet_bTagWeightLFStats1Up,   "Jet_bTagWeightLFStats1Up/F");
  globalTree->Branch("Jet_bTagWeightLFStats1Down", &Jet_bTagWeightLFStats1Down, "Jet_bTagWeightLFStats1Down/F");
  globalTree->Branch("Jet_bTagWeightLFStats2Up",   &Jet_bTagWeightLFStats2Up,   "Jet_bTagWeightLFStats2Up/F");
  globalTree->Branch("Jet_bTagWeightLFStats2Down", &Jet_bTagWeightLFStats2Down, "Jet_bTagWeightLFStats2Down/F");

  globalTree->Branch("Jet_bTagWeightHFUp",         &Jet_bTagWeightHFUp,        "Jet_bTagWeightHFUp/F");
  globalTree->Branch("Jet_bTagWeightHFDown",       &Jet_bTagWeightHFDown,      "Jet_bTagWeightHFDown/F");
  globalTree->Branch("Jet_bTagWeightHFStats1Up",   &Jet_bTagWeightHFStats1Up,   "Jet_bTagWeightHFStats1Up/F");
  globalTree->Branch("Jet_bTagWeightHFStats1Down", &Jet_bTagWeightHFStats1Down, "Jet_bTagWeightHFStats1Down/F");
  globalTree->Branch("Jet_bTagWeightHFStats2Up",   &Jet_bTagWeightHFStats2Up,   "Jet_bTagWeightHFStats2Up/F");
  globalTree->Branch("Jet_bTagWeightHFStats2Down", &Jet_bTagWeightHFStats2Down, "Jet_bTagWeightHFStats2Down/F");

  globalTree->Branch("Jet_bTagWeightcErr1Up",    &Jet_bTagWeightcErr1Up,   "Jet_bTagWeightcErr1Up/F");
  globalTree->Branch("Jet_bTagWeightcErr1Down",  &Jet_bTagWeightcErr1Down, "Jet_bTagWeightcErr1Down/F");
  globalTree->Branch("Jet_bTagWeightcErr2Up",    &Jet_bTagWeightcErr2Up,   "Jet_bTagWeightcErr2Up/F");
  globalTree->Branch("Jet_bTagWeightcErr2Down",  &Jet_bTagWeightcErr2Down, "Jet_bTagWeightcErr2Down/F");

  globalTree->Branch("Jet_btagCSV",&Jet_btagCSV, "Jet_btagCSV/F"); //Pending patch
  globalTree->Branch("Jet_btagCMVA",&Jet_btagCMVA, "Jet_btagCMVA/F");
  globalTree->Branch("Jet_btagCSVV0",&Jet_btagCSVV0, "Jet_btagCSVV0/F");
  globalTree->Branch("Jet_btagCMVAV2",&Jet_btagCMVAV2, "Jet_btagCMVAV2/F");

  globalTree->Branch("Jet_mcPt",&Jet_mcPt, "Jet_mcPt/F");

  globalTree->Branch("Jet_corr",         &Jet_corr, "Jet_corr/F");
  globalTree->Branch("Jet_corr_JECUp",   &Jet_corr_JECUp, "Jet_corr_JECUp/F");
  globalTree->Branch("Jet_corr_JECDown", &Jet_corr_JECDown, "Jet_corr_JECDown/F");

  globalTree->Branch("Jet_corr_JER",     &Jet_corr_JER, "Jet_corr_JER");
  globalTree->Branch("Jet_corr_JERUp",   &Jet_corr_JERUp, "Jet_corr_JERUp/F");
  globalTree->Branch("Jet_corr_JERDown", &Jet_corr_JERDown, "Jet_corr_JERDown/F");
 
  globalTree->Branch("nPVs", &_nPv, "_nPv/I");    

  if(!runningOnData_){
   // globalTree->Branch("Genb_pt",&Genb_pt);
   // globalTree->Branch("Genb_eta",&Genb_eta);
   // globalTree->Branch("Genb_phi",&Genb_phi);
   // globalTree->Branch("Genb_mass",&Genb_mass);

    globalTree->Branch("Genb_pt",   &Genb_pt,   "Genb_pt/F");
    globalTree->Branch("Genb_eta",  &Genb_eta,  "Genb_eta/F");
    globalTree->Branch("Genb_phi",  &Genb_phi,  "Genb_phi/F");
    globalTree->Branch("Genb_mass", &Genb_mass, "Genb_mass/F");

    globalTree->Branch("nTrueInt",      &npT, "npT/F");
    globalTree->Branch("InTime_PileUp", &npIT, "npIT/F");   //Not yet included with proper var/branch name

    globalTree->Branch("Jet_mcFlavour",&Jet_mcFlavour, "Jet_mcFlavour/I");
  }

  globalTree->Branch("Jet_partonFlavour", &Jet_partonFlavour, "Jet_partonFlavour/I");
  globalTree->Branch("Jet_hadronFlavour", &Jet_hadronFlavour, "Jet_hadronFlavour/I");

  globalTree->Branch("nGenJet",&nGenJet, "nGenJet/I");                    //Pending patch
  globalTree->Branch("Jet_mcMatchId",&Jet_mcMatchId, "Jet_mcMatchId/F");

  globalTree->Branch("HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v",&HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v, "HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v/I");
  globalTree->Branch("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v",&HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v, "HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v/I");
  globalTree->Branch("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v",&HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v, "HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v/I");
  globalTree->Branch("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v",&HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v, "HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v/I");

  globalTree->Branch("HLT_HH4bAll", &HLT_HH4bAll, "HLT_HH4bAll/I");

  globalTree->Branch("genWeight",   &genWeight, "genWeight/F");

  globalTree->Branch("puWeight",     &puWeight,     "puWeight/F");
  globalTree->Branch("puWeightUp",   &puWeightUp,   "puWeightUp/F");
  globalTree->Branch("puWeightDown", &puWeightDown, "puWeightDown/F");
  
  globalTree->Branch("xsec",        &xsec,        "xsec/F");
  globalTree->Branch("json",        &json,        "json/F");
  globalTree->Branch("json_silver", &json_silver, "json_silver/F");

  globalTree->Branch("bTagWeight_LFUp",&bTagWeight_LFUp, "bTagWeight_LFUp/F");
  globalTree->Branch("bTagWeight_LFDown",&bTagWeight_LFDown, "bTagWeight_LFDown/F");
  globalTree->Branch("bTagWeight_LFStats1Up",&bTagWeight_LFStats1Up, "bTagWeight_LFStats1Up/F");
  globalTree->Branch("bTagWeight_LFStats1Down",&bTagWeight_LFStats1Down, "bTagWeight_LFStats1Down/F");
  globalTree->Branch("bTagWeight_LFStats2Down",&bTagWeight_LFStats2Down, "bTagWeight_LFStats2Down/F");
  globalTree->Branch("bTagWeight_LFStats2Up",&bTagWeight_LFStats2Up, "bTagWeight_LFStats2Up/F");

  globalTree->Branch("bTagWeight_HFUp",&bTagWeight_HFUp, "bTagWeight_HFUp/F");
  globalTree->Branch("bTagWeight_HFDown",&bTagWeight_HFDown, "bTagWeight_HFDown/F");
  globalTree->Branch("bTagWeight_HFStats1Up",&bTagWeight_HFStats1Up, "bTagWeight_HFStats1Up/F");
  globalTree->Branch("bTagWeight_HFStats1Down",&bTagWeight_HFStats1Down, "bTagWeight_HFStats1Down/F");
  globalTree->Branch("bTagWeight_HFStats2Up",&bTagWeight_HFStats2Up, "bTagWeight_HFStats2Up/F");
  globalTree->Branch("bTagWeight_HFStats2Down",&bTagWeight_HFStats2Down, "bTagWeight_HFStats2Down/F");

  globalTree->Branch("bTagWeight_cErr1Down",&bTagWeight_cErr1Down, "bTagWeight_cErr1Down/F");
  globalTree->Branch("bTagWeight_cErr1Up",&bTagWeight_cErr1Up, "bTagWeight_cErr1Up/F");
  globalTree->Branch("bTagWeight_cErr2Up",&bTagWeight_cErr2Up, "bTagWeight_cErr2Up/F");
  globalTree->Branch("bTagWeight_cErr2Down",&bTagWeight_cErr2Down, "bTagWeight_cErr2Down/F");

  globalTree->Branch("bTagWeight_JESDown",&bTagWeight_JESDown, "bTagWeight_JESDown/F");
  globalTree->Branch("bTagWeight_JESUp",&bTagWeight_JESUp, "bTagWeight_JESUp/F");
 
  globalTree->Branch("met_pt",   &met_pt,   "met_pt/F");
  globalTree->Branch("met_eta",  &met_eta,  "met_eta/F");
  globalTree->Branch("met_phi",  &met_phi,  "met_phi/F");
  globalTree->Branch("met_mass", &met_mass, "met_mass/F");

  globalTree->Branch("GenJet_pt",   &GenJet_pt);
  globalTree->Branch("GenJet_eta",  &GenJet_eta);
  globalTree->Branch("GenJet_phi",  &GenJet_phi);
  globalTree->Branch("GenJet_mass", &GenJet_mass);

  globalTree->Branch("LHE_weights_scale_wgt",&LHE_weights_scale_wgt, "LHE_weights_scale_wgt/F");
  globalTree->Branch("nLHE_weights_pdf",&nLHE_weights_pdf, "nLHE_weights_pdf/F");
  globalTree->Branch("LHE_weights_pdf_id",&LHE_weights_pdf_id, "LHE_weights_pdf_id/F");
  globalTree->Branch("LHE_weights_pdf_wgt",&LHE_weights_pdf_wgt, "LHE_weights_pdf_wgt/F");

}

//define this as a plug-in
DEFINE_FWK_MODULE(HAA4bAnalysis);
