import ROOT
import math
import numpy as np

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H500_A200")

##Global constants
#PT1_MIN = 165.
#PT2_MIN = 130.
#PT3_MIN = 130.
#PT4_MIN = 110.
#DELTA_PHI_MIN = 0.9
#DELTA_ETA_MAX = 1.5
#PT_PAIR1_MIN = 150.
#PT_PAIR2_MIN = 150.

PT1_MIN = 165.
PT2_MIN = 130.
PT3_MIN = 130.
PT4_MIN = 110.
DELTA_PHI_MIN = 0.9
DELTA_ETA_MAX = 1.5
PT_PAIR1_MIN = 0.
PT_PAIR2_MIN = 0.

##Normalize to this luminsity, in fb-1
luminosity_norm = 10.
##Make signal histos larger
signal_magnify = 10.

list_histos = ["h_jet1pt", "h_jet2pt", "h_jet3pt", "h_jet4pt", "h_delta_Phi_pair", "h_delta_Eta_pair", "h_pt_pair1", "h_pt_pair2", "h_m4b"]

def select_all_but_one(cutstring):

    selection_bools = dict()
    selection_bools["h_jet1pt"] = jet1pt > PT1_MIN
    selection_bools["h_jet2pt"] = jet2pt > PT2_MIN
    selection_bools["h_jet3pt"] = jet3pt > PT3_MIN
    selection_bools["h_jet4pt"] = jet4pt > PT4_MIN
    selection_bools["h_delta_Phi_pair"] = delta_Phi > DELTA_PHI_MIN
    selection_bools["h_delta_Eta_pair"] = abs(delta_Eta) < DELTA_ETA_MAX
    selection_bools["h_pt_pair1"] = pt_pair1 > PT_PAIR1_MIN
    selection_bools["h_pt_pair2"] = pt_pair2 > PT_PAIR1_MIN

    result = True

    for hname in selection_bools:
        if cutstring == hname:
            continue
        else:
            result = result and selection_bools[hname]

    return result

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

##Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##Get the handlers for all the histos
hs = dict()
h_QCD = dict()
h_base = dict()
h_Signal = dict()

for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")

##Define the histos to be created
for sample_name in samplename_list:
    if "QCD" in sample_name:
        continue

    h_base[sample_name+list_histos[0]]  = ROOT.TH1F(sample_name+list_histos[0], "p_{t} of 1st jet", 50, PT1_MIN, 500.)
    h_base[sample_name+list_histos[1]]  = ROOT.TH1F(sample_name+list_histos[1], "p_{t} of 2nd jet", 50, PT2_MIN, 500.)
    h_base[sample_name+list_histos[2]]  = ROOT.TH1F(sample_name+list_histos[2], "p_{t} of 3rd jet", 50, PT3_MIN, 500.)
    h_base[sample_name+list_histos[3]]  = ROOT.TH1F(sample_name+list_histos[3], "p_{t} of 4th jet", 50, PT4_MIN, 500.)
    h_base[sample_name+list_histos[4]]  = ROOT.TH1F(sample_name+list_histos[4], "#Delta_{#phi} of the two jet pairs", 30, 0., 6.28)
    h_base[sample_name+list_histos[5]]  = ROOT.TH1F(sample_name+list_histos[5], "#Delta_{#eta} of the two jet pairs", 20, -5., 5.)
    h_base[sample_name+list_histos[6]]  = ROOT.TH1F(sample_name+list_histos[6], "p_T of the first jet pair", 50, 0., 500.)
    h_base[sample_name+list_histos[7]]  = ROOT.TH1F(sample_name+list_histos[7], "p_T of the second jet pair", 50, 0., 500.)
    h_base[sample_name+list_histos[8]]  = ROOT.TH1F(sample_name+list_histos[8], "4jets invariant mass", 100, 0., 1000.)

h_QCD[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{t} of 1st jet", 50, PT1_MIN, 500.)
h_QCD[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{t} of 2nd jet", 50, PT2_MIN, 500.)
h_QCD[list_histos[2]]  = ROOT.TH1F(list_histos[2], "p_{t} of 3rd jet", 50, PT3_MIN, 500.)
h_QCD[list_histos[3]]  = ROOT.TH1F(list_histos[3], "p_{t} of 4th jet", 50, PT4_MIN, 500.)
h_QCD[list_histos[4]]  = ROOT.TH1F(list_histos[4], "#Delta_{#phi} of the two jet pairs", 30, 0., 6.28)
h_QCD[list_histos[5]]  = ROOT.TH1F(list_histos[5], "#Delta_{#eta} of the two jet pairs", 20, -5., 5.)
h_QCD[list_histos[6]]  = ROOT.TH1F(list_histos[6], "p_T of the first jet pair", 50, 0., 500.)
h_QCD[list_histos[7]]  = ROOT.TH1F(list_histos[7], "p_T of the second jet pair", 50, 0., 500.)
h_QCD[list_histos[8]]  = ROOT.TH1F(list_histos[8], "4jets invariant mass", 100, 0., 1000.)

##Graphics stuff
canvas = dict()
for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname)
leg1 = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg1.SetHeader("Samples considered")

##Loop on samples, and then on events, and merge QCD stuff
idx_sample = 0
for name_sample in samplename_list:

    if not "QCD" in name_sample:
        QCDflag = False
    else:
        QCDflag = True

    norm_factor = Norm_Map[name_sample]*luminosity_norm
    mytree = root_file[name_sample].Get("HZZ4bAnalysis/mytree")

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break

        nb = mytree.GetEntry( jentry )
        if nb <= 0:
            continue

        jet1_4mom = mytree.jet1_4mom
        jet2_4mom = mytree.jet2_4mom
        jet3_4mom = mytree.jet3_4mom
        jet4_4mom = mytree.jet4_4mom

        jet1pt = jet1_4mom.Pt()
        jet2pt = jet2_4mom.Pt()
        jet3pt = jet3_4mom.Pt()
        jet4pt = jet4_4mom.Pt()

        combination_flag = myWF.get_best_combination(jet1_4mom,jet2_4mom,jet3_4mom,jet4_4mom);

        if combination_flag == 1:
            p_pair1 = jet1_4mom + jet2_4mom
            p_pair2 = jet3_4mom + jet4_4mom
        elif combination_flag == 2:
            p_pair1 = jet1_4mom + jet3_4mom
            p_pair2 = jet2_4mom + jet4_4mom
        elif combination_flag == 3:
            p_pair1 = jet1_4mom + jet4_4mom
            p_pair2 = jet2_4mom + jet3_4mom

        delta_Phi = abs(p_pair1.Phi() - p_pair2.Phi())
        delta_Eta = p_pair1.Eta() - p_pair2.Eta()

        pt_pair1 = p_pair1.Pt()
        pt_pair2 = p_pair2.Pt()

        totaljets_4mom = p_pair1 + p_pair2
        m4b = totaljets_4mom.M()

        if QCDflag:
            if select_all_but_one("h_jet1pt"):
                h_QCD["h_jet1pt"].Fill(jet1pt,norm_factor)
            if select_all_but_one("h_jet2pt"):
                h_QCD["h_jet2pt"].Fill(jet2pt,norm_factor)
            if select_all_but_one("h_jet3pt"):
                h_QCD["h_jet3pt"].Fill(jet3pt,norm_factor)
            if select_all_but_one("h_jet4pt"):
                h_QCD["h_jet4pt"].Fill(jet4pt,norm_factor)
            if select_all_but_one("h_delta_Phi_pair"):
                h_QCD["h_delta_Phi_pair"].Fill(delta_Phi,norm_factor)
            if select_all_but_one("h_delta_Eta_pair"):
                h_QCD["h_delta_Eta_pair"].Fill(delta_Eta,norm_factor)
            if select_all_but_one("h_pt_pair1"):
                h_QCD["h_pt_pair1"].Fill(pt_pair1,norm_factor)
            if select_all_but_one("h_pt_pair2"):
                h_QCD["h_pt_pair2"].Fill(pt_pair2,norm_factor)
            if select_all_but_one("h_m4b"):
                h_QCD["h_m4b"].Fill(m4b,norm_factor)
        else:
            if select_all_but_one("h_jet1pt"):
                h_base[name_sample+"h_jet1pt"].Fill(jet1pt,norm_factor)
            if select_all_but_one("h_jet2pt"):
                h_base[name_sample+"h_jet2pt"].Fill(jet2pt,norm_factor)
            if select_all_but_one("h_jet3pt"):
                h_base[name_sample+"h_jet3pt"].Fill(jet3pt,norm_factor)
            if select_all_but_one("h_jet4pt"):
                h_base[name_sample+"h_jet4pt"].Fill(jet4pt,norm_factor)
            if select_all_but_one("h_delta_Phi_pair"):
                h_base[name_sample+"h_delta_Phi_pair"].Fill(delta_Phi,norm_factor)
            if select_all_but_one("h_delta_Eta_pair"):
                h_base[name_sample+"h_delta_Eta_pair"].Fill(delta_Eta,norm_factor)
            if select_all_but_one("h_pt_pair1"):
                h_base[name_sample+"h_pt_pair1"].Fill(pt_pair1,norm_factor)
            if select_all_but_one("h_pt_pair2"):
                h_base[name_sample+"h_pt_pair2"].Fill(pt_pair2,norm_factor)
            if select_all_but_one("h_m4b"):
                h_base[name_sample+"h_m4b"].Fill(m4b,norm_factor)

    if QCDflag:
        continue

    for idx_histo,hname in enumerate(list_histos):
        h_base[name_sample+hname].SetFillColor(idx_sample+3)
        if name_sample == myWF.sig_samplename:
                h_base[name_sample+hname].SetLineStyle(2)   #dashed
                h_base[name_sample+hname].SetLineColor(4)   #blue
                h_base[name_sample+hname].SetLineWidth(4)   #kind of thik
        if idx_histo == 0:
            leg1.AddEntry(h_base[name_sample+hname],name_sample,"f")

        hs[hname].Add(h_base[name_sample+hname])

    idx_sample += 1
    if idx_sample == 1:
        idx_sample += 1

for idx_histo,hname in enumerate(list_histos):
    h_QCD[hname].SetFillColor(2) #Red
    if idx_histo == 0:
        leg1.AddEntry(h_QCD[hname],"QCD","f")
    hs[hname].Add(h_QCD[hname])


for hname in list_histos:
    canvas[hname].cd()
    hs[hname].Draw("histo")
    h_base[myWF.sig_samplename+hname].Scale(signal_magnify)
    h_base[myWF.sig_samplename+hname].Draw("SAMEHISTO")
    leg1.Draw()
    canvas[hname].SaveAs("plots/tree_" + hname + ".gif")
