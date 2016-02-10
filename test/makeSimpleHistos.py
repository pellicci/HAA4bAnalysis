import ROOT

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H500_A200")

##Do all the scaling to this luminosity, in fb-1
luminosity_norm = 10.
signal_magnify = 100.

##These are the histograms to be merged
list_histos = ["h_jet1pt", "h_jet2pt", "h_jet3pt", "h_jet4pt", "h_jet1Btag", "h_jet2Btag", "h_jet3Btag", "h_jet4Btag", "h_delta_Phi_pair", "h_delta_Eta_pair", "h_Events", "h_jet1eta", "h_jet2eta"]

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

##Get the handlers for all the histos
h_QCD = dict()
h_QCD[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{t} of 1st jet", 200, 0., 500.)
h_QCD[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{t} of 2nd jet", 200, 0., 500.)
h_QCD[list_histos[2]]  = ROOT.TH1F(list_histos[2], "p_{t} of 3rd jet", 200, 0., 500.)
h_QCD[list_histos[3]]  = ROOT.TH1F(list_histos[3], "p_{t} of 4th jet", 200, 0., 500.)
h_QCD[list_histos[4]]  = ROOT.TH1F(list_histos[4], "Btag value of 1st jet", 50, 0.89, 1.)
h_QCD[list_histos[5]]  = ROOT.TH1F(list_histos[5], "Btag value of 2nd jet", 50, 0.89, 1.)
h_QCD[list_histos[6]]  = ROOT.TH1F(list_histos[6], "Btag value of 3rd jet", 50, 0.89, 1.)
h_QCD[list_histos[7]]  = ROOT.TH1F(list_histos[7], "Btag value of 4th jet", 50, 0.89, 1.)
h_QCD[list_histos[8]]  = ROOT.TH1F(list_histos[8], "#Delta_{#phi} between the two jet pairs", 50, 0., 6.28)
h_QCD[list_histos[9]]  = ROOT.TH1F(list_histos[9], "#Delta_{#eta} between the two jet pairs", 50, -10., 10.)
h_QCD[list_histos[10]] = ROOT.TH1F(list_histos[10], "Events after pre-selection steps", 6, 0., 6.)
h_QCD[list_histos[11]] = ROOT.TH1F(list_histos[11], "#eta of 1st jet", 50, -10., 10.)
h_QCD[list_histos[12]] = ROOT.TH1F(list_histos[12], "#eta of 2nd jet", 50, -10., 10.)

hs = dict()
h_signal = dict()
for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")
    h_signal[hname] = h_QCD[hname].Clone()

##Graphics stuff
canvas = dict()
for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname)
leg1 = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg1.SetHeader("Samples considered")

##collect all the root files
samplename_list = myWF.get_samples_names(False)
root_file = myWF.get_root_files(False)

##Loop on all the histos to be created
for idx_hname,hname in enumerate(list_histos): 

    ##Merge the QCD samples
    for name_sample in samplename_list:
        if not "QCD" in name_sample:
            continue

        norm_factor = Norm_Map[name_sample]*luminosity_norm

        h_tmp = ROOT.TH1F(root_file[name_sample].Get("HZZ4bAnalysis/" + hname))
        h_tmp.Scale(norm_factor)

        h_QCD[hname].Add(h_tmp)

    h_QCD[hname].SetFillColor(2) #Red
    if idx_hname == 0:
        leg1.AddEntry(h_QCD[hname],"QCD","f")

    hs[hname].Add(h_QCD[hname])

    ##Loop on all the other background samples
    idx_color = 0
    for name_sample in samplename_list:
        if "QCD" in name_sample:
            continue

        norm_factor = Norm_Map[name_sample]*luminosity_norm
        h_tmp = ROOT.TH1F(root_file[name_sample].Get("HZZ4bAnalysis/" + hname))

        h_tmp.Scale(norm_factor)
        h_tmp.SetFillColor(idx_color+3)
        idx_color += 1
        if idx_hname == 0:
            leg1.AddEntry(h_tmp,name_sample,"f")

        hs[hname].Add(h_tmp)

##Now deal with the signal
root_file_signal = ROOT.TFile(myWF.sig_filename)
for idx_hname,hname in enumerate(list_histos): 
    norm_factor = Norm_Map[myWF.sig_samplename]*luminosity_norm
    h_tmp = ROOT.TH1F(root_file_signal.Get("HZZ4bAnalysis/" + hname))

    h_tmp.Scale(norm_factor)
    h_tmp.SetLineStyle(2)   #dashed
    h_tmp.SetLineColor(4)   #blue
    h_tmp.SetLineWidth(4)   #kind of thik
    if idx_hname == 0:
        leg1.AddEntry(h_tmp,myWF.sig_samplename,"f")

    h_signal[hname] = h_tmp.Clone(hname)

    hs[hname].Add(h_tmp)

##Put lables in the nEvents plot

##Manual tuning of plots


for hname in list_histos:
    canvas[hname].cd()
    hs[hname].Draw("histo")
    if hname == "h_Events" or hname == "h_delta_Phi_pair":
        canvas[hname].SetLogy(1)

    h_signal[hname].Scale(signal_magnify)
    h_signal[hname].Draw("same")
    leg1.Draw()
    canvas[hname].SaveAs("plots/" + hname + ".gif")
