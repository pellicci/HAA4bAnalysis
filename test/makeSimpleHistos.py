import ROOT
import os

norm_filename = "Normalizations_table.txt"
dir_back_input = "rootfiles/backgrounds/"

dir_sig_input = "rootfiles/signals/"
sig_samplename = "Signal_H500_A200"
sig_filename = dir_sig_input + "HAA4bAnalysis_" + sig_samplename + ".root"

##Do all the scaling to this luminosity, in fb-1
luminosity_norm = 10.
signal_magnify = 100000.

##These are the histograms to be merged
list_histos = ["h_jet1pt", "h_jet2pt", "h_jet3pt", "h_jet4pt", "h_jet1Btag", "h_jet2Btag", "h_jet3Btag", "h_jet4Btag", "h_delta_Phi_pair", "h_delta_Eta_pair"]

def get_normalizations_map(filename):
    in_file = open(filename,"r")
    norm_map = dict()

    for line in in_file:
        data_norm = line.split()
        norm_map[data_norm[0]] = float(data_norm[1])

    return norm_map


##Here starts the program
Norm_Map = get_normalizations_map(norm_filename)

##Get the handlers for all the histos
hs = dict()
for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")

h_QCD = dict()
h_QCD[list_histos[0]] = ROOT.TH1F(list_histos[0], "p_{t} of 1st jet", 200, 0., 500.)
h_QCD[list_histos[1]] = ROOT.TH1F(list_histos[1], "p_{t} of 2nd jet", 200, 0., 500.)
h_QCD[list_histos[2]] = ROOT.TH1F(list_histos[2], "p_{t} of 3rd jet", 200, 0., 500.)
h_QCD[list_histos[3]] = ROOT.TH1F(list_histos[3], "p_{t} of 4th jet", 200, 0., 500.)
h_QCD[list_histos[4]] = ROOT.TH1F(list_histos[4], "Btag value of 1st jet", 50, 0., 1.)
h_QCD[list_histos[5]] = ROOT.TH1F(list_histos[5], "Btag value of 2nd jet", 50, 0., 1.)
h_QCD[list_histos[6]] = ROOT.TH1F(list_histos[6], "Btag value of 3rd jet", 50, 0., 1.)
h_QCD[list_histos[7]] = ROOT.TH1F(list_histos[7], "Btag value of 4th jet", 50, 0., 1.)
h_QCD[list_histos[8]] = ROOT.TH1F(list_histos[8], "#Delta_{#phi} between the two jet pairs", 30, 0., 3.14)
h_QCD[list_histos[9]] = ROOT.TH1F(list_histos[9], "#Delta_{#eta} between the two jet pairs", 50, -10., 10.)

h_signal = dict()
h_signal[list_histos[0]] = ROOT.TH1F(list_histos[0] + "_Signal", "p_{t} of 1st jet", 200, 0., 500.)
h_signal[list_histos[1]] = ROOT.TH1F(list_histos[1] + "_Signal", "p_{t} of 2nd jet", 200, 0., 500.)
h_signal[list_histos[2]] = ROOT.TH1F(list_histos[2] + "_Signal", "p_{t} of 3rd jet", 200, 0., 500.)
h_signal[list_histos[3]] = ROOT.TH1F(list_histos[3] + "_Signal", "p_{t} of 4th jet", 200, 0., 500.)
h_signal[list_histos[4]] = ROOT.TH1F(list_histos[4] + "_Signal", "Btag value  of 1st jet", 50, 0., 1.)
h_signal[list_histos[5]] = ROOT.TH1F(list_histos[5] + "_Signal", "Btag value  of 2nd jet", 50, 0., 1.)
h_signal[list_histos[6]] = ROOT.TH1F(list_histos[6] + "_Signal", "Btag value  of 3rd jet", 50, 0., 1.)
h_signal[list_histos[7]] = ROOT.TH1F(list_histos[7] + "_Signal", "Btag value  of 4th jet", 50, 0., 1.)
h_signal[list_histos[8]] = ROOT.TH1F(list_histos[8] + "_Signal", "#Delta_{#phi} between the two jet pairs", 30, 0., 3.14)
h_signal[list_histos[9]] = ROOT.TH1F(list_histos[9] + "_Signal", "#Delta_{#eta} between the two jet pairs", 50, -10., 10.)


##Graphics stuff
canvas = dict()
for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname)
leg1 = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg1.SetHeader("Samples considered")

##collect all the root files
list_dirs_bkg = os.listdir(dir_back_input)
root_file = dict()
samplename_list = []
for dirname in list_dirs_bkg:

    tmp_samplename = dirname.split("HAA4bAnalysis_")[1]
    tmp_samplename2 = tmp_samplename.replace(".root","")
    samplename_list.append(tmp_samplename2)
    root_file[tmp_samplename2] = ROOT.TFile(dir_back_input + dirname)

##Merge the QCD samples
for idx_hname,hname in enumerate(list_histos): 

    for name_sample in samplename_list:
        if not "QCD" in name_sample:
            continue

        norm_factor = Norm_Map[name_sample]*luminosity_norm

        h_tmp = ROOT.TH1F(root_file[name_sample].Get("HZZ4bAnalysis/" + hname))
        h_tmp.Scale(norm_factor)
        h_QCD[hname].Add(h_tmp,norm_factor)

    h_QCD[hname].SetFillColor(2) #Red
    if idx_hname == 0:
        leg1.AddEntry(h_QCD[hname],"QCD","f")

    hs[hname].Add(h_QCD[hname])

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
root_file_signal = ROOT.TFile(sig_filename)
for idx_hname,hname in enumerate(list_histos): 
    norm_factor = Norm_Map[sig_samplename]*luminosity_norm
    h_tmp = ROOT.TH1F(root_file_signal.Get("HZZ4bAnalysis/" + hname))
    h_tmp.Scale(norm_factor)

    h_tmp.SetLineStyle(2)   #dashed
    h_tmp.SetLineColor(4)   #blue
    h_tmp.SetLineWidth(4)   #kind of thik
    if idx_hname == 0:
        leg1.AddEntry(h_tmp,sig_samplename,"f")

    h_signal[hname] = h_tmp.Clone(hname)
    hs[hname].Add(h_tmp)

for hname in list_histos:
    canvas[hname].cd()
    hs[hname].Draw()
    h_signal[hname].Scale(signal_magnify)
    h_signal[hname].Draw("same")
    leg1.Draw()
    canvas[hname].SaveAs("plots/" + hname + ".gif")
