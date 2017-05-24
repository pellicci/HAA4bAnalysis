
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

  // vectors to store MC information for background analysis
  std::vector<float> Genb_pt;
  std::vector<float> Genb_phi;
  std::vector<float> Genb_eta;
  std::vector<float> Genb_mass;

  //Met variables
  float gMet_pt;
  float gMet_phi;
  float gMet_eta;
  float gMet_mass;
  float Jet_btag;

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
  double bjet_scaleFactor;
  double bjet_scaleFactor_ScaleUp;
  double bjet_scaleFactor_ScaleDown;

  //Heavy Flavour (HF) Scale factors
  float HF_scaleFactor;
  float HF_scaleFactor_ScaleUp;
  float HF_scaleFactor_ScaleDown;

  //Light Flavor Scale (LF) factors
  float LF_scaleFactor;
  float LF_scaleFactor_ScaleUp;
  float LF_scaleFactor_ScaleDown;
 
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


  float xsec;

  float puWeight;
  float puWeightUp; 
  float puWeightDown; 

  float genWeight;
  float json;
  float json_silver;

  float bTagWeight_LFUp;
  float bTagWeight_LFDown;

  float bTagWeight_LFStats1Up;
  float bTagWeight_LFStats1Down;
  float bTagWeight_LFStats2Down;
  float bTagWeight_LFStats2Up;

  float bTagWeight_HFUp;
  float bTagWeight_HFDown;

  float bTagWeight_HFStats1Up;
  float bTagWeight_HFStats1Down;
  float bTagWeight_HFStats2Up;
  float bTagWeight_HFStats2Down;

  float bTagWeight_cErr1Down;
  float bTagWeight_cErr1Up;
  float bTagWeight_cErr2Up;
  float bTagWeight_cErr2Down;

  float bTagWeight_JESDown;
  float bTagWeight_JESUp;

  float bTagWeight;

  float Jet_bTagWeight;

  float Jet_bTagWeightJESUp;
  float Jet_bTagWeightJESDown;

  float Jet_bTagWeightLFUp;   //start
  float Jet_bTagWeightLFDown;

  float Jet_bTagWeightLFStats1Up;
  float Jet_bTagWeightLFStats1Down;
  float Jet_bTagWeightLFStats2Up;
  float Jet_bTagWeightLFStats2Down;

  float Jet_bTagWeightHFUp;
  float Jet_bTagWeightHFDown;

  float Jet_bTagWeightHFStats1Up;
  float Jet_bTagWeightHFStats1Down;
  float Jet_bTagWeightHFStats2Up;
  float Jet_bTagWeightHFStats2Down;

  float Jet_bTagWeightcErr1Up;
  float Jet_bTagWeightcErr1Down;
  float Jet_bTagWeightcErr2Up;
  float Jet_bTagWeightcErr2Down; 

  float Jet_btagCSV;
  float Jet_btagCMVA; 
  float Jet_btagCSVV0; 
  float Jet_btagCMVAV2; 

  float Jet_mcPt;

  float Jet_corr; 
  float Jet_corr_JECUp; 
  float Jet_corr_JECDown; 

  float Jet_corr_JER;
  float Jet_corr_JERUp; 
  float Jet_corr_JERDown; 

  float met_pt;
  float met_eta; 
  float met_phi; 
  float met_mass; 

  float Jet_pt;
  float Jet_eta; 
  float Jet_phi; 
  float Jet_mass; 

  float GenJet_pt; 
  float GenJet_eta; 
  float GenJet_phi; 
  float GenJet_mass; 

  float LHE_weights_scale_wgt; 
  float nLHE_weights_pdf;
  float LHE_weights_pdf_id; 
  float LHE_weights_pdf_wgt;

 //Btagging Info of Jets 
  float var_jet1Btag;
  float var_jet2Btag;
  float var_jet3Btag;
  float var_jet4Btag;

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
