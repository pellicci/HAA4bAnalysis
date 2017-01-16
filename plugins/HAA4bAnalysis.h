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

  const edm::InputTag jets_;
  const edm::InputTag globaljets_; 
  const edm::InputTag genParticles_;
  const edm::InputTag met_;
  const edm::InputTag globalmet_;
  std::string bdiscr_;
  double minPt_high_;
  double minPt_low_;
  double minCSV_;
  bool runningOnData_;
  const edm::InputTag pvCollection_;  
  const edm::InputTag bsCollection_;  
  const edm::InputTag PileupSrc_;

  edm::LumiReWeighting Lumiweights_; 
  edm::Service<TFileService> fs;

  void create_Histos_and_Trees();
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
  TH1F* h_nPv;      // No. of primary verticies histogram
  TH1F* h_Rw_nPv;   //Reweighted No. of primary verticies

  TH1F* h_PUInTime;  //In time PileUp
  TH1F* h_PUTrue;    //True number of PileUp
  TH1F* h_PUWeight;         //Weight   
  TH1F* h_Rw_PUTrue;     //Reweighted True number of PileUp Interactions
  TH1F* h_Rw_PUInTime;   //Reweighted Intime PileUp
  
  //TTree stuff
  TTree *mytree; //this is for the "general" analysis
  TTree *globalTree; //this is for the study of QCD background

  TLorentzVector *jet1_4mom_tree;
  TLorentzVector *jet2_4mom_tree;
  TLorentzVector *jet3_4mom_tree;
  TLorentzVector *jet4_4mom_tree;
  TLorentzVector *jet1_4mom_tree_fit;
  TLorentzVector *jet2_4mom_tree_fit;
  TLorentzVector *jet3_4mom_tree_fit;
  TLorentzVector *jet4_4mom_tree_fit;

  //Vectors to store Jet Information
  std::vector<float> Jet_pt;
  std::vector<float> Jet_phi;
  std::vector<float> Jet_eta;
  std::vector<float> Jet_mass;
  std::vector<float> Jet_btag;
  std::vector<float> jetbTagWeight;
  std::vector<int>   Jet_hadflavrs;
  std::vector<int>   Jet_partnflavrs;
  std::vector<float> Jet_corr;
  std::vector<float> Jet_corr_shifted;

  // vectors to store MC information for background analysis
  std::vector<float> Genb_pt;
  std::vector<float> Genb_phi;
  std::vector<float> Genb_eta;
  std::vector<float> Genb_mass;
  std::vector<int>   Genb_flavor;
 
 //Hadron and parton Flavour Information vectors
  std::vector<int> Genb_hadflavrs;
  std::vector<int> Genb_partnflavrs;

  // vectors to store global MET
  std::vector<float> gMet_pt;
  std::vector<float> gMet_phi;
  std::vector<float> gMet_eta;
  std::vector<float> gMet_mass;
 // std::vector<float> Jet_btag;

  //Event Info like event and run number, lumi sec., 
  //std::vector<edm::EventID> event_nmbr;
  std::vector<int> event_nmbr;
  std::vector<int> run_nmbr;
  std::vector<int> lumi_blck;
  std::vector<bool> is_data;
  std::vector<bool> is_json;
  std::vector<bool> is_json_silver;
  std::vector<float> is_xsec;
  std::vector<float> pu_weight;
  std::vector<float> pu_weightUp;
  std::vector<float> pu_weightDown;

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
    

  unsigned int _Noutputs; //to limit the couts of debug info

  //Tokens
  edm::EDGetTokenT<std::vector<pat::Jet> > jetstoken_; 
  edm::EDGetTokenT<std::vector<pat::Jet> > globaljetstoken_; 
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlestoken_; 
  edm::EDGetTokenT<std::vector<pat::MET> > metToken_;
  edm::EDGetTokenT<std::vector<pat::MET> > globalmetToken_;
  edm::EDGetTokenT<reco::VertexCollection> tok_Vertex_; 
  edm::EDGetTokenT<reco::BeamSpot>         tok_beamspot_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupSummaryToken_;
};
