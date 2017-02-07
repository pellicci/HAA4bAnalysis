import ROOT

fIn = ROOT.TFile("inputs/Nominal_training.root")
tree_Signal     = fIn.Get("tree_signal")
tree_Background = fIn.Get("tree_background")

fOut = ROOT.TFile("outputs/Nominal_training.root","RECREATE")

ROOT.TMVA.Tools.Instance()

factory = ROOT.TMVA.Factory("TMVAClassification", fOut,":".join(["!V","!Silent","Color","DrawProgressBar","Transformations=I;D;P;G,D","AnalysisType=Classification"]))

factory.AddVariable("j1_pt","F")
factory.AddVariable("j2_pt","F")
factory.AddVariable("j3_pt","F")
factory.AddVariable("j4_pt","F")

factory.AddSpectator("j1_btag","F")
factory.AddSpectator("j2_btag","F")
factory.AddSpectator("j3_btag","F")
factory.AddSpectator("j4_btag","F")

factory.AddVariable("delta_phi","F")
factory.AddVariable("delta_eta","F")

factory.AddVariable("abs_massRatio_jetpair","F")

factory.AddSignalTree    (tree_Signal)
factory.AddBackgroundTree(tree_Background)

factory.SetSignalWeightExpression("evt_weight")
factory.SetBackgroundWeightExpression("evt_weight");

mycuts = ROOT.TCut("j1_btag > 0.8 && j2_btag > 0.8 && j3_btag > 0.8 && j4_btag > 0.8")
mycutb = ROOT.TCut("j1_btag > 0.8 && j2_btag > 0.8 && j3_btag > 0.8 && j4_btag > 0.8")

factory.PrepareTrainingAndTestTree(mycuts, mycutb, ":".join(["SplitMode=random","!V"]) )

method_cuts = factory.BookMethod(ROOT.TMVA.Types.kCuts,"Cuts",":".join(["!H","!V","FitMethod=MC","EffSel","SampleSize=10000000","VarProp=FSmart"]))

#method_btd  = factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", ":".join(["!H","!V","NTrees=1000", "MinNodeSize=2.5%","MaxDepth=3","BoostType=AdaBoost","AdaBoostBeta=0.5","SeparationType=GiniIndex","nCuts=20"]))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()

fOut.Close()
