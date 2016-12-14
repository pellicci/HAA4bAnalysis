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
//
//


// system include files
#include <memory>
#include <algorithm>

//ROOT includes
#include <TH1F.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TTree.h>

//CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// Kinematic fitter
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"

// vertex inclusions
#include "DataFormats/VertexReco/interface/Vertex.h" 
#include "DataFormats/BeamSpot/interface/BeamSpot.h" 
// pileUp inclusions
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// Jet Energy corrections
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/JetCorrector/interface/JetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

typedef math::XYZTLorentzVector LorentzVector;

#include "HAA4bAnalysis.h"

// constructors and destructor
HAA4bAnalysis::HAA4bAnalysis(const edm::ParameterSet& iConfig) :
  jets_(iConfig.getParameter<edm::InputTag>("jets")), 
  globaljets_(iConfig.getParameter<edm::InputTag>("globaljets")), 
  genParticles_(iConfig.getParameter<edm::InputTag>("genParticles")), 
  bdiscr_(iConfig.getParameter<std::string>("BTagAlgo")),
  minPt_high_(iConfig.getParameter<double>("minPt_high")),
  minPt_low_(iConfig.getParameter<double>("minPt_low")),
  minCSV_(iConfig.getParameter<double>("minCSV")),
  runningOnData_(iConfig.getParameter<bool>("runningOnData")),   
  pvCollection_(iConfig.getParameter<edm::InputTag>("pvCollection")),   
  bsCollection_(iConfig.getParameter<edm::InputTag>("bsCollection")),  
  PileupSrc_(iConfig.getParameter<edm::InputTag>("PileupSrc")) //,
{
  jetstoken_          = consumes<std::vector<pat::Jet> >(jets_); 
  globaljetstoken_    = consumes<std::vector<pat::Jet> >(globaljets_);
  genParticlestoken_  = consumes<std::vector<reco::GenParticle> >(genParticles_);
  tok_Vertex_         = consumes<reco::VertexCollection>(pvCollection_);  
  tok_beamspot_       = consumes<reco::BeamSpot>(edm::InputTag(bsCollection_));
  pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag(PileupSrc_)); 

   //Few Counters
  _Nevents_processed = 0.;
  _Nevents_4jets     = 0.;
  _Nevents_4bjets    = 0.;
  _Nevents_ptpass    = 0.;
  _Nevents_mpairs    = 0.;
  _Nevents_deltaM    = 0.;
  _Nevents_passed    = 0.;

  _Noutputs = 0;

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
  // no point to continue if there aren't 4 jets
  if(jets->size() < 4) return;

  // get the Handle of the primary vertex collection and remove the beamspot
  edm::Handle<reco::BeamSpot> bmspot;
  iEvent.getByLabel(bsCollection_,bmspot);

  edm::Handle<reco::VertexCollection> pvcoll;  
  iEvent.getByLabel(pvCollection_,pvcoll);

  // require in the event that there is at least one reconstructed vertex
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
  PU_Weight = 1.;

  if (!runningOnData_){
    std::cout<<"Running on Mote Carlo "<<std::endl;
    edm::Handle<std::vector< PileupSummaryInfo>>  PupInfo;
    iEvent.getByLabel(PileupSrc_, PupInfo);  
  
    std::vector<PileupSummaryInfo>::const_iterator PVI; 
 
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      int BX = PVI->getBunchCrossing();
      if(BX == 0) {
	npT = PVI->getTrueNumInteractions();
	npIT = PVI->getPU_NumInteractions();
      }
    } 

    // calculate weight using above code
    PU_Weight = Lumiweights_.weight(npT);
    if(_show_output) std::cout<<"PU_Weight for MC is "<<PU_Weight<<std::endl;
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
  // define a jet handle and get the jets without special selections
  edm::Handle<std::vector<pat::Jet> > globaljets;
  iEvent.getByLabel(globaljets_, globaljets);

  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  if(!runningOnData_) iEvent.getByLabel(genParticles_, genParticles);
  fill_global_Tree(globaljets, genParticles);

  //Select the highest csv jets
  pat::Jet jet1 = jets->at(nj1);
  pat::Jet jet2 = jets->at(nj2);
  pat::Jet jet3 = jets->at(nj3);
  pat::Jet jet4 = jets->at(nj4);

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
  if (!runningOnData_) Lumiweights_=edm::LumiReWeighting("MC_Recent_25ns_2015.root", "pileUpData_fromJson.root", "pileup", "pileup");
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

  std::cout << "#######################" << std::endl;
  std::cout << "EVENT SUMMARY" << std::endl;
  std::cout << "Number of events processed = " << _Nevents_processed << std::endl;
  std::cout << "Number of events with 4 jets = " << _Nevents_4jets << std::endl;
  std::cout << "Number of events with 3 bjets = " << _Nevents_4bjets << std::endl;
  std::cout << "Number of events passing pt cuts = " << _Nevents_ptpass << std::endl;
  std::cout << "Number of events with m_12/34 > 120 GeV = " << _Nevents_mpairs << std::endl;
  std::cout << "Number of events with delta_m_12/34 < 100 GeV = " << _Nevents_deltaM << std::endl;
  std::cout << "Number of events passing the selection = " << _Nevents_passed << std::endl;
  std::cout << "#######################" << std::endl;
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
  fitter.addMeasParticle( &jet1 );
  fitter.addMeasParticle( &jet2 );
  fitter.addMeasParticle( &jet3 );
  fitter.addMeasParticle( &jet4 );
  fitter.addConstraint( &mCons1 );

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
    std::cout << "=============================================" << std ::endl;
    std::cout << "-> Number of measured Particles  : " << fitter.nbMeasParticles() << std::endl;
    std::cout << "-> Number of unmeasured particles: " << fitter.nbUnmeasParticles() << std::endl;
    std::cout << "-> Number of constraints         : " << fitter.nbConstraints() << std::endl;
    std::cout << "-> Number of degrees of freedom  : " << fitter.getNDF() << std::endl;
    std::cout << "-> Number of parameters A        : " << fitter.getNParA() << std::endl;
    std::cout << "-> Number of parameters B        : " << fitter.getNParB() << std::endl;
    std::cout << "-> Maximum number of iterations  : " << fitter.getMaxNumberIter() << std::endl;
    std::cout << "-> Maximum deltaS                : " << fitter.getMaxDeltaS() << std::endl;
    std::cout << "-> Maximum F                     : " << fitter.getMaxF() << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std ::endl;
    std::cout << "-> Status                        : " << fitter.getStatus() << std::endl;
    std::cout << "-> Number of iterations          : " << fitter.getNbIter() << std::endl;
    std::cout << "-> S                             : " << fitter.getS() << std::endl;
    std::cout << "-> F                             : " << fitter.getF() << std::endl;
    std::cout << "=============================================" << std ::endl;
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

void HAA4bAnalysis::fill_global_Tree(edm::Handle<std::vector<pat::Jet> >& globaljets, edm::Handle<std::vector<reco::GenParticle> >& genParticles){

  //Clean up the vectors to contain the jet information
  Jet_pt.clear();
  Jet_eta.clear();
  Jet_phi.clear();
  Jet_mass.clear();
  Jet_btag.clear();

  for(auto jet = globaljets->begin(); jet != globaljets->end(); ++jet){
    float thept = jet->p4().Pt();
    float theeta = jet->p4().Eta();
    float thephi = jet->p4().Phi();
    float themass = jet->p4().M();
    float thecsv = jet->bDiscriminator(bdiscr_);

    Jet_pt.push_back(thept);
    Jet_eta.push_back(theeta);
    Jet_phi.push_back(thephi);
    Jet_mass.push_back(themass);
    Jet_btag.push_back(thecsv);
  }

  Genb_pt.clear();
  Genb_eta.clear();
  Genb_phi.clear();
  Genb_mass.clear();

  if(!runningOnData_){
    //std::cout << "Number of genparticles = " << genParticles->size() << std::endl;

    for(auto genpart = genParticles->begin(); genpart != genParticles->end(); genpart++){

      //std::cout << "Id of genparticle = " << genpart->pdgId() << std::endl;
      if( abs(genpart->pdgId()) != 5) continue;
      if( abs(genpart->mother(0)->pdgId()) != 36) continue; //A -> bb

      Genb_pt.push_back(genpart->pt());
      Genb_eta.push_back(genpart->eta());
      Genb_phi.push_back(genpart->phi());
      Genb_mass.push_back(genpart->mass());
    }
  }

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

  h_nPv = fs->make<TH1F>("h_nPv", "No. of primary verticies", 50, 0., 50.);              //Few more Histogram inclusions
  h_Rw_nPv = fs->make<TH1F>("h_Rw_nPv", "Reweighted No. of primary verticies", 50, 0., 50.);
  h_PUInTime= fs->make<TH1F>("h_PUInTime","Input No. in-time pileup interactions",50,0.,50.);
  h_PUTrue= fs->make<TH1F>("h_PUTrue","Input True pileup interactions",50,0.,50.);
  h_Rw_PUTrue = fs->make<TH1F>("h_Rw_PUTrue","Reweighted True pileup interactions",50,0.,50.);
  h_Rw_PUInTime = fs->make<TH1F>("h_Rw_PUInTime","Reweighted in-time pileup interactions",50,0.,50.);
  h_PUWeight = fs->make<TH1F>("h_PUWeight","Event weight",50,0.,50.);
 //h_WeightVsNint = fs->make<TProfile>("h_WeightVsNint","Event weight vs N_int",50,0.,50.,0.,10.);

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

  // pileUp histograms
  mytree->Branch("PUTrue", &npT, "npT/F");
  mytree->Branch("PUInTime", &npIT, "npIT/F");
  mytree->Branch("PUWeight", &PU_Weight, "PU_Weight/F");

  globalTree = fs->make<TTree>("globalTree", "Tree containing most of the event info to study background");
  globalTree->Branch("Jet_pt",&Jet_pt);
  globalTree->Branch("Jet_eta",&Jet_eta);
  globalTree->Branch("Jet_phi",&Jet_phi);
  globalTree->Branch("Jet_mass",&Jet_mass);
  globalTree->Branch("Jet_btag",&Jet_btag);
  if(!runningOnData_){
    globalTree->Branch("Genb_pt",&Genb_pt);
    globalTree->Branch("Genb_eta",&Genb_eta);
    globalTree->Branch("Genb_phi",&Genb_phi);
    globalTree->Branch("Genb_mass",&Genb_mass);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(HAA4bAnalysis);
