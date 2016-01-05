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

//CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

typedef math::XYZTLorentzVector LorentzVector;

//
// class declaration
//

class HAA4bAnalysis : public edm::EDAnalyzer {
   public:
      explicit HAA4bAnalysis(const edm::ParameterSet&);
      ~HAA4bAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      int get_best_combination(LorentzVector m1, LorentzVector m2, LorentzVector m3, LorentzVector m4);
      bool check_combinations(LorentzVector m1, LorentzVector m2, LorentzVector m3, LorentzVector m4);

      // ----------member data ---------------------------
      const edm::InputTag jets_;
      std::string bdiscr_;
      double minPt_;

      edm::Service<TFileService> fs;

      TH1F* h_jet1pt;
      TH1F* h_jet2pt;
      TH1F* h_jet3pt;
      TH1F* h_jet4pt;

      TH1F* h_jet1Btag;
      TH1F* h_jet2Btag;
      TH1F* h_jet3Btag;
      TH1F* h_jet4Btag;

      TH1F* h_m_pair1;
      TH1F* h_m_pair2;
      TH1F* h_m_4b;

      TH1F* h_delta_Phi_pair;
      TH1F* h_delta_Eta_pair;

      int _Nevents_processed;
      int _Nevents_passed;
  
};


// constructors and destructor
HAA4bAnalysis::HAA4bAnalysis(const edm::ParameterSet& iConfig) :
  jets_(iConfig.getParameter<edm::InputTag>("jets")),
  bdiscr_(iConfig.getParameter<std::string>("BTagAlgo")),
  minPt_(iConfig.getParameter<double>("minPt"))
{
  _Nevents_processed = 0;
  _Nevents_passed = 0;

  h_jet1pt = fs->make<TH1F>("h_jet1pt", "P_t of the 1st jet", 200, 0., 500.);
  h_jet2pt = fs->make<TH1F>("h_jet2pt", "P_t of the 2nd jet", 200, 0., 500.);
  h_jet3pt = fs->make<TH1F>("h_jet3pt", "P_t of the 3rd jet", 200, 0., 500.);
  h_jet4pt = fs->make<TH1F>("h_jet4pt", "P_t of the 4th jet", 200, 0., 500.);

  h_jet1Btag = fs->make<TH1F>("h_jet1Btag", "B-tag discriminant of the 1st jet", 200, 0., 500.);
  h_jet2Btag = fs->make<TH1F>("h_jet2Btag", "B-tag discriminant of the 2nd jet", 200, 0., 500.);
  h_jet3Btag = fs->make<TH1F>("h_jet3Btag", "B-tag discriminant of the 3rd jet", 200, 0., 500.);
  h_jet4Btag = fs->make<TH1F>("h_jet4Btag", "B-tag discriminant of the 4th jet", 200, 0., 500.);

  h_m_pair1 = fs->make<TH1F>("h_m_pair1", "Invariant mass of the first jet pair", 400, 120., 800.);
  h_m_pair2 = fs->make<TH1F>("h_m_pair2", "Invariant mass of the second jet pair", 400, 120., 800.);
  h_m_4b    = fs->make<TH1F>("h_m_4b", "Invariant mass of the four b jets", 400, 120., 1000.);

  h_delta_Phi_pair = fs->make<TH1F>("h_delta_Phi_pair", "#Delta_{#phi} between the two jet pairs", 30, -3.14, 3.14);
  h_delta_Eta_pair = fs->make<TH1F>("h_delta_Eta_pair", "#Delta_{#eta} between the two jet pairs", 50, -10., 10.);

}


HAA4bAnalysis::~HAA4bAnalysis()
{
}



// member functions

// ------------ method called for each event  ------------
void HAA4bAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  _Nevents_processed++;

  // define a jet handle
  edm::Handle<std::vector<pat::Jet> > jets;
  // get jets from the event
  iEvent.getByLabel(jets_, jets);

  if(jets->size() < 4) return;

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
  for(auto jet = jets->begin(); jet != jets->end(); ++jet){
      float thept = jet->p4().Pt();

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

  pat::Jet jet1 = jets->at(nj1);
  pat::Jet jet2 = jets->at(nj2);
  pat::Jet jet3 = jets->at(nj3);
  pat::Jet jet4 = jets->at(nj4);

  LorentzVector jet1_4mom = jet1.p4();
  LorentzVector jet2_4mom = jet2.p4();
  LorentzVector jet3_4mom = jet3.p4();
  LorentzVector jet4_4mom = jet4.p4();
 
  if(jet4_4mom.Pt() < minPt_) return;

  //Refuse to continue if no combination is above 120 GeV
  if( !check_combinations(jet1_4mom,jet2_4mom,jet3_4mom,jet4_4mom) ) return;

  //Convention:
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

  //Avoid a case in which the difference between the masses is above 100 GeV
  if(fabs (p_pair1.M() - p_pair2.M() ) > 100.) return;

  //Angular distributions between the two pairs
  float delta_Phi_pairs = p_pair1.Phi() - p_pair2.Phi();
  float delta_Eta_pairs = p_pair1.Eta() - p_pair2.Eta();

  float total_Mass_4b = (p_pair1 + p_pair2).M();
  if(total_Mass_4b < 250.) return;

  //Now plot a few quantities for the jets
  h_jet1pt->Fill(jet1.p4().Pt());
  h_jet2pt->Fill(jet2.p4().Pt());
  h_jet3pt->Fill(jet3.p4().Pt());
  h_jet4pt->Fill(jet4.p4().Pt());

  h_jet1Btag->Fill(jet1.bDiscriminator(bdiscr_));
  h_jet2Btag->Fill(jet2.bDiscriminator(bdiscr_));
  h_jet3Btag->Fill(jet3.bDiscriminator(bdiscr_));
  h_jet4Btag->Fill(jet4.bDiscriminator(bdiscr_));

  //And the invariant masses
  h_m_pair1->Fill(p_pair1.M());
  h_m_pair2->Fill(p_pair2.M());
  h_m_4b->Fill(total_Mass_4b);

  //And the angular "distributions"
  h_delta_Phi_pair->Fill(delta_Phi_pairs);
  h_delta_Eta_pair->Fill(delta_Eta_pairs);

}


// ------------ method called once each job just before starting event loop  ------------
void HAA4bAnalysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void HAA4bAnalysis::endJob() 
{
  std::cout << "EVENT SUMMARY" << std::endl;
  std::cout << "Number of events processed = " << _Nevents_processed << std::endl;
  
}

int HAA4bAnalysis::get_best_combination(LorentzVector m1, LorentzVector m2, LorentzVector m3, LorentzVector m4){

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

bool HAA4bAnalysis::check_combinations(LorentzVector m1, LorentzVector m2, LorentzVector m3, LorentzVector m4){

  bool b12 = (m1+m2).M() > 120.;
  bool b13 = (m1+m3).M() > 120.;
  bool b14 = (m1+m4).M() > 120.;
  bool b34 = (m3+m4).M() > 120.;
  bool b24 = (m2+m4).M() > 120.;
  bool b23 = (m2+m3).M() > 120.;

  return b12 || b13 || b14 || b34 || b24 || b23;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HAA4bAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HAA4bAnalysis);
