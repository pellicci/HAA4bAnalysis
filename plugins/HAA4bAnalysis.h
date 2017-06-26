
//---------- class declaration----------

class HAA4bAnalysis : public edm::EDAnalyzer {
public:
  explicit HAA4bAnalysis(const edm::ParameterSet&);
  ~HAA4bAnalysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  const edm::InputTag   jets_;
  const edm::InputTag   globaljets_; 
  const edm::InputTag   genParticles_;
  const edm::InputTag   met_;
  const edm::InputTag   globalmet_;
        std::string     bdiscr_;
  double                minPt_high_;
  double                minPt_low_;
  double                minCSV_;
  bool                  runningOnData_;
  const  edm::InputTag  pvCollection_;  
  const  edm::InputTag  bsCollection_;  
  const  edm::InputTag  PileupSrc_;
  const  edm::InputTag  rhoSrc_;

  //For PileUp Stuff: Nominal, Up and Down
  edm::LumiReWeighting Lumiweights_;
  edm::LumiReWeighting LumiweightsUp_;
  edm::LumiReWeighting LumiweightsDown_;

  std::vector<std::string> jecPayloadNames_;
  std::string              jecUncName_;
  boost::shared_ptr<JetCorrectionUncertainty> jecUnc_;
  boost::shared_ptr<FactorizedJetCorrector>   jec_;
 
  //File Service 
  edm::Service<TFileService> fs;

  void create_Histos_and_Trees();

 //tree function to fill the background varaiables
  void fill_global_Tree(edm::Handle<std::vector<pat::Jet> >& globaljets, edm::Handle<std::vector<reco::GenParticle> >& genParticles, edm::Handle<std::vector<pat::MET> > &globalmets);

  int get_best_combination(LorentzVector& m1, LorentzVector& m2, LorentzVector& m3, LorentzVector& m4);
  bool check_combinations(LorentzVector& m1, LorentzVector& m2, LorentzVector& m3, LorentzVector& m4, float mcut);
  TKinFitter get_fitted_candidate(pat::Jet& Jet1, pat::Jet& Jet2, pat::Jet& Jet3, pat::Jet& Jet4, int best_combination);

  // ----------member data ---------------------------
  TH1F* h_jet1pt;
  TH1F* h_jet2pt;
  TH1F* h_jet3pt;
  TH1F* h_jet4pt;

  TH1F* h_jet1eta;
  TH1F* h_jet2eta;
  TH1F* h_jet3eta;
  TH1F* h_jet4eta;

  TH1F* h_jet1Btag;
  TH1F* h_jet2Btag;
  TH1F* h_jet3Btag;
  TH1F* h_jet4Btag;

  TH1F* h_m_pair1;
  TH1F* h_m_pair2;
  TH1F* h_m_4b;
  TH1F* h_m_4b_fitted;
  TH2F* h_m4b_m12;
  TH2F* h_m4b_m34;

  TH1F* h_delta_Phi_pair;
  TH1F* h_delta_Eta_pair;

  TH1F* h_Events;

  TH1F* h_nPv;          // No. of primary verticies histogram
  TH1F* h_Rw_nPv;       //Reweighted No. of primary verticies
  TH1F* h_PUInTime;     //In time PileUp
  TH1F* h_PUTrue;       //True number of PileUp
  TH1F* h_PUWeight;     //Weight   
  TH1F* h_Rw_PUTrue;    //Reweighted True number of PileUp Interactions
  TH1F* h_Rw_PUInTime;  //Reweighted Intime PileUp
  
  //TTree stuff
  TTree *mytree;        //this is for the "general" analysis
  TTree *globalTree;    //this is for the study of QCD background

  TLorentzVector *jet1_4mom_tree;
  TLorentzVector *jet2_4mom_tree;
  TLorentzVector *jet3_4mom_tree;
  TLorentzVector *jet4_4mom_tree;
  TLorentzVector *jet1_4mom_tree_fit;
  TLorentzVector *jet2_4mom_tree_fit;
  TLorentzVector *jet3_4mom_tree_fit;
  TLorentzVector *jet4_4mom_tree_fit;
//==============================================

   const int Max_Jets = 200;

  // vectors to store MC information for background analysis
  std::vector<float> Genb_pt;
  std::vector<float> Genb_phi;
  std::vector<float> Genb_eta;
  std::vector<float> Genb_mass;

  //Met variables
 /* float gMet_pt;
  float gMet_phi;
  float gMet_eta;
  float gMet_mass;
  //float Jet_btag;
*/
  //Jet Resolutions variables
  float Jet_resolution;
  float Jet_resolutionSF;
  float Jet_resolutionSF_Up;
  float Jet_resolutionSF_Down;

  JME::JetParameters            Parameters;
  JME::JetResolution            Resolution;
  JME::JetResolutionScaleFactor SF_Resolution;

 //JEC uncertainity variables
  double corr;       //JEC factor
  double corrUp;     //JEC ScaleUp
  double corrDown;   //JEC ScaleDown 
  double JEC_Uncertainity;
  double ptCor_ScaleUp;
  double ptCor_ScaleDown;

  //Btag Weight and Scale variation
  float* bjet_scaleFactor           = new float[Max_Jets]();
  float* bjet_scaleFactor_ScaleUp   = new float[Max_Jets]();
  float* bjet_scaleFactor_ScaleDown = new float[Max_Jets]();

  //Event Info like event and run number, lumi sec., 
  ULong64_t evt;
  unsigned int run;
  unsigned int lumi;

  Int_t isData;
  Int_t Jet_id;
  Int_t Jet_mcFlavour;
  Int_t Jet_partonFlavour;
  Int_t Jet_hadronFlavour;
  Int_t nGenJet;
  Int_t Jet_mcMatchId;

  Int_t HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v;
  Int_t HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v;
  Int_t HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v;
  Int_t HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v;
  Int_t HLT_HH4bAll;
  
 //int Jet_Size_it;

   int Jet_iterator;
   int globalJet_iterator;
   int genJet_iterator;
   //int genbJet_iterator;
   int met_iterator;

  //std::vector<int> Jet_Size_V;

  float xsec;

  float puWeight;
  float puWeightUp; 
  float puWeightDown; 

  float genWeight;
  float json;
  float json_silver;

/*  float* bTagWeight_LFUp   = new float[Max_Jets];
  float* bTagWeight_LFDown = new float[Max_Jets];

  float* bTagWeight_LFStats1Up   = new float[Max_Jets];
  float* bTagWeight_LFStats1Down = new float[Max_Jets];
  float* bTagWeight_LFStats2Down = new float[Max_Jets];
  float* bTagWeight_LFStats2Up   = new float[Max_Jets];

  float* bTagWeight_HFUp   = new float[Max_Jets];
  float* bTagWeight_HFDown = new float[Max_Jets];

  float* bTagWeight_HFStats1Up   = new float[Max_Jets];
  float* bTagWeight_HFStats1Down = new float[Max_Jets];
  float* bTagWeight_HFStats2Up   = new float[Max_Jets];
  float* bTagWeight_HFStats2Down = new float[Max_Jets];

  float* bTagWeight_cErr1Down = new float[Max_Jets];
  float* bTagWeight_cErr1Up   = new float[Max_Jets];
  float* bTagWeight_cErr2Up   = new float[Max_Jets];
  float* bTagWeight_cErr2Down = new float[Max_Jets];

  float* bTagWeight_JESDown = new float[Max_Jets] ;
  float* bTagWeight_JESUp   = new float[Max_Jets];

  float* bTagWeight = new float[Max_Jets];
*/

  float* Jet_bTagWeight = new float[Max_Jets]();

  float* Jet_bTagWeightJESUp = new float[Max_Jets]();
  float* Jet_bTagWeightJESDown  = new float[Max_Jets]();

  float* Jet_bTagWeightLFUp   = new float[Max_Jets]();   //start
  float* Jet_bTagWeightLFDown = new float[Max_Jets]();

  float* Jet_bTagWeightLFStats1Up   = new float[Max_Jets]();
  float* Jet_bTagWeightLFStats1Down = new float[Max_Jets]();
  float* Jet_bTagWeightLFStats2Up   = new float[Max_Jets]();
  float* Jet_bTagWeightLFStats2Down = new float[Max_Jets]();

  float* Jet_bTagWeightHFUp   = new float[Max_Jets]();
  float* Jet_bTagWeightHFDown = new float[Max_Jets]();

  float* Jet_bTagWeightHFStats1Up   = new float[Max_Jets]();
  float* Jet_bTagWeightHFStats1Down = new float[Max_Jets]();
  float* Jet_bTagWeightHFStats2Up   = new float[Max_Jets]();
  float* Jet_bTagWeightHFStats2Down = new float[Max_Jets]();

  float* Jet_bTagWeightcErr1Up   = new float[Max_Jets]();
  float* Jet_bTagWeightcErr1Down = new float[Max_Jets]();
  float* Jet_bTagWeightcErr2Up   = new float[Max_Jets]();
  float* Jet_bTagWeightcErr2Down = new float[Max_Jets](); 


  float* Jet_btagCSV    = new float[Max_Jets]();
  float* Jet_btagCMVA   = new float[Max_Jets](); 
  float* Jet_btagCSVV0  = new float[Max_Jets](); 
  float* Jet_btagCMVAV2 = new float[Max_Jets](); 

  float* Jet_mcPt = new float[Max_Jets]();

  float* Jet_corr         = new float[Max_Jets](); 
  float* Jet_corr_JECUp   = new float[Max_Jets](); 
  float* Jet_corr_JECDown = new float[Max_Jets](); 

  float* Jet_corr_JER     = new float[Max_Jets]();
  float* Jet_corr_JERUp   = new float[Max_Jets](); 
  float* Jet_corr_JERDown = new float[Max_Jets](); 

  float* met_pt   = new float[Max_Jets]();
  float* met_eta  = new float[Max_Jets](); 
  float* met_phi  = new float[Max_Jets](); 
  float* met_mass = new float[Max_Jets](); 

 
  float* Jet_pt   = new float[Max_Jets]();
  float* Jet_eta  = new float[Max_Jets](); 
  float* Jet_phi  = new float[Max_Jets](); 
  float* Jet_mass = new float[Max_Jets](); 
  float* Jet_btag = new float[Max_Jets]();

  float* GenJet_pt    = new float[Max_Jets](); 
  float* GenJet_eta   = new float[Max_Jets](); 
  float* GenJet_phi   = new float[Max_Jets](); 
  float* GenJet_mass  = new float[Max_Jets](); 

  float LHE_weights_scale_wgt; 
  float nLHE_weights_pdf;
  float LHE_weights_pdf_id; 
  float LHE_weights_pdf_wgt;

 //Btagging Info of Jets 
  float var_jet1Btag;
  float var_jet2Btag;
  float var_jet3Btag;
  float var_jet4Btag;

// Input for Gamma Variables
  float jet1_energy;
  float jet2_energy;
  float jet3_energy;
  float jet4_energy;

  float jet1_mass;
  float jet2_mass;
  float jet3_mass;
  float jet4_mass;


  //A few counters
  float _Nevents_processed;
  float _Nevents_4jets;
  float _Nevents_4bjets;
  float _Nevents_ptpass;
  float _Nevents_mpairs;
  float _Nevents_deltaM;
  float _Nevents_passed;
  unsigned int _nPv; 

  //Pile  Up Info used in "mytree"
  float npT;
  float npIT;
  float PU_Weight;   
  float PU_WeightUp;
  float PU_WeightDown;  
  // Few variables for background estimation and initialization
  int   is_data;
  int   is_not_data;
  float is_mc;
  float is_not_mc;
  float is_json;
  float is_silver_json;
  float is_xsec;
  float is_trigger; 
  int   is_mcFlavour;
  int   nGenJets;
  int   jet_id;
  int   jet_mcMatchId;

  unsigned int is_run;
  unsigned int _Noutputs; //to limit the couts of debug info

  //Tokens
  edm::EDGetTokenT<std::vector<pat::Jet> >            jetstoken_; 
  edm::EDGetTokenT<std::vector<pat::Jet> >            globaljetstoken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> >   genParticlestoken_; 
  edm::EDGetTokenT<std::vector<pat::MET> >            metToken_;
  edm::EDGetTokenT<std::vector<pat::MET> >            globalmetToken_;
  edm::EDGetTokenT<reco::VertexCollection>            tok_Vertex_; 
  edm::EDGetTokenT<reco::BeamSpot>                    tok_beamspot_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>>    pileupSummaryToken_;
  edm::EDGetTokenT<double>                            rhoSrc_token_;
};
