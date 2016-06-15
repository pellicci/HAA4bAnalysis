import ROOT
import os
import math
import numpy as np

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H800_A300")

##Global constants
##This is for loose preselection
PT1_MIN = 110.     #100
PT2_MIN = 70.     #100
PT3_MIN = 70.     #70
PT4_MIN = 60.
ETA1_MAX = 2.5
ETA2_MAX = 2.5
ETA3_MAX = 2.5
ETA4_MAX = 2.5
DELTA_PHI_MIN = 0.2
PT_PAIR1_MIN = 20.
PT_PAIR2_MIN = 20.
M_DIJET_MIN = 130.
JET1_BTAG = 0.89 #0.97
JET2_BTAG = 0.89 #0.97
JET3_BTAG = 0.89 #0.97
JET4_BTAG = 0.89 #0.97

##Normalize to this luminsity, in fb-1
luminosity_norm = 5.
##Make signal histos larger
signal_magnify = 50.

output_dir = "plots"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

list_histos = ["h_jet1pt", "h_jet2pt", "h_jet3pt", "h_jet4pt", "h_delta_Phi_pair", "h_delta_Eta_pair", "h_pt_pair1", "h_pt_pair2", "h_m4b", "h_jet1eta", "h_jet2eta", "h_jet3eta", "h_jet4eta", "h_jet1Btag", "h_jet2Btag", "h_jet3Btag", "h_jet4Btag", "h_m12", "h_m34","h_nPv" ]

def select_all_but_one(cutstring):

    selection_bools = dict()
    selection_bools["h_jet1pt"] = jet1pt > PT1_MIN
    selection_bools["h_jet2pt"] = jet2pt > PT2_MIN
    selection_bools["h_jet3pt"] = jet3pt > PT3_MIN
    selection_bools["h_jet4pt"] = jet4pt > PT4_MIN

    selection_bools["h_jet1eta"] = jet1eta < ETA1_MAX
    selection_bools["h_jet2eta"] = jet2eta < ETA2_MAX
    selection_bools["h_jet3eta"] = jet3eta < ETA3_MAX
    selection_bools["h_jet4eta"] = jet4eta < ETA4_MAX

    selection_bools["h_delta_Phi_pair"] = delta_Phi > DELTA_PHI_MIN
    selection_bools["h_delta_Eta_pair"] = abs(delta_Eta) < DELTA_ETA_MAX
    selection_bools["h_pt_pair1"] = pt_pair1 > PT_PAIR1_MIN
    selection_bools["h_pt_pair2"] = pt_pair2 > PT_PAIR1_MIN

    selection_bools["h_jet1Btag"] = mytree.jet1Btag > JET1_BTAG
    selection_bools["h_jet2Btag"] = mytree.jet2Btag > JET2_BTAG
    selection_bools["h_jet3Btag"] = mytree.jet3Btag > JET3_BTAG
    selection_bools["h_jet4Btag"] = mytree.jet4Btag > JET4_BTAG

    selection_bools["h_m12"] = m_pair1 > M_DIJET_MIN
    selection_bools["h_m34"] = m_pair2 > M_DIJET_MIN

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
    if "Signal" in sample_name and not sample_name == myWF.sig_samplename:
        continue

    h_base[sample_name+list_histos[0]]  = ROOT.TH1F(sample_name+list_histos[0], "p_{t} of 1st jet", 25, 50., 500.)
    h_base[sample_name+list_histos[1]]  = ROOT.TH1F(sample_name+list_histos[1], "p_{t} of 2nd jet", 25, 50., 500.)
    h_base[sample_name+list_histos[2]]  = ROOT.TH1F(sample_name+list_histos[2], "p_{t} of 3rd jet", 25, 50., 500.)
    h_base[sample_name+list_histos[3]]  = ROOT.TH1F(sample_name+list_histos[3], "p_{t} of 4th jet", 25, 50., 500.)
    h_base[sample_name+list_histos[4]]  = ROOT.TH1F(sample_name+list_histos[4], "#Delta_{#phi} of the two jet pairs", 15, 0., 3.14)
    h_base[sample_name+list_histos[5]]  = ROOT.TH1F(sample_name+list_histos[5], "#Delta_{#eta} of the two jet pairs", 20, -5., 5.)
    h_base[sample_name+list_histos[6]]  = ROOT.TH1F(sample_name+list_histos[6], "p_T of the first jet pair", 50, 0., 500.)
    h_base[sample_name+list_histos[7]]  = ROOT.TH1F(sample_name+list_histos[7], "p_T of the second jet pair", 50, 0., 500.)
    h_base[sample_name+list_histos[8]]  = ROOT.TH1F(sample_name+list_histos[8], "4jets invariant mass", 40, 0., 1000.)
    h_base[sample_name+list_histos[9]]  = ROOT.TH1F(sample_name+list_histos[9], "#eta of 1st jet", 20, -7., 7.)
    h_base[sample_name+list_histos[10]] = ROOT.TH1F(sample_name+list_histos[10], "#eta of 2nd jet", 20, -7., 7.)
    h_base[sample_name+list_histos[11]] = ROOT.TH1F(sample_name+list_histos[11], "#eta of 3rd jet", 20, -7., 7.)
    h_base[sample_name+list_histos[12]] = ROOT.TH1F(sample_name+list_histos[12], "#eta of 4th jet", 20, -7., 7.)

    h_base[sample_name+list_histos[13]] = ROOT.TH1F(sample_name+list_histos[13], "B-tag of 1st jet", 13, 0.6, 1.)
    h_base[sample_name+list_histos[14]] = ROOT.TH1F(sample_name+list_histos[14], "B-tag of 2nd jet", 13, 0.6, 1.)
    h_base[sample_name+list_histos[15]] = ROOT.TH1F(sample_name+list_histos[15], "B-tag of 3rd jet", 13, 0.6, 1.)
    h_base[sample_name+list_histos[16]] = ROOT.TH1F(sample_name+list_histos[16], "B-tag of 4th jet", 13, 0., 1.)

    h_base[sample_name+list_histos[17]] = ROOT.TH1F(sample_name+list_histos[17], "m_{12} invariant mass", 15, 0., 600.)
    h_base[sample_name+list_histos[18]] = ROOT.TH1F(sample_name+list_histos[18], "m_{34} invariant mass", 15, 0., 600.)
    h_base[sample_name+list_histos[19]] = ROOT.TH1F(sample_name+list_histos[19], "No. of primary verticies", 50, 0., 50.)

h_QCD[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{t} of 1st jet", 25, 50., 500.)
h_QCD[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{t} of 2nd jet", 25, 50., 500.)
h_QCD[list_histos[2]]  = ROOT.TH1F(list_histos[2], "p_{t} of 3rd jet", 25, 50., 500.)
h_QCD[list_histos[3]]  = ROOT.TH1F(list_histos[3], "p_{t} of 4th jet", 25, 50., 500.)
h_QCD[list_histos[4]]  = ROOT.TH1F(list_histos[4], "#Delta_{#phi} of the two jet pairs", 15, 0., 3.14)
h_QCD[list_histos[5]]  = ROOT.TH1F(list_histos[5], "#Delta_{#eta} of the two jet pairs", 20, -5., 5.)
h_QCD[list_histos[6]]  = ROOT.TH1F(list_histos[6], "p_T of the first jet pair", 50, 0., 500.)
h_QCD[list_histos[7]]  = ROOT.TH1F(list_histos[7], "p_T of the second jet pair", 50, 0., 500.)
h_QCD[list_histos[8]]  = ROOT.TH1F(list_histos[8], "4jets invariant mass", 40, 0., 1000.)
h_QCD[list_histos[9]]  = ROOT.TH1F(list_histos[9], "#eta of 1st jet", 20, -7., 7.)
h_QCD[list_histos[10]] = ROOT.TH1F(list_histos[10], "#eta of 2nd jet", 20, -7., 7.)
h_QCD[list_histos[11]] = ROOT.TH1F(list_histos[11], "#eta of 3rd jet", 20, -7., 7.)
h_QCD[list_histos[12]] = ROOT.TH1F(list_histos[12], "#eta of 4th jet", 20, -7., 7.)

h_QCD[list_histos[13]] = ROOT.TH1F(list_histos[13], "B-tag of 1st jet", 13, 0.6, 1.)
h_QCD[list_histos[14]] = ROOT.TH1F(list_histos[14], "B-tag of 2nd jet", 13, 0.6, 1.)
h_QCD[list_histos[15]] = ROOT.TH1F(list_histos[15], "B-tag of 3rd jet", 13, 0.6, 1.)
h_QCD[list_histos[16]] = ROOT.TH1F(list_histos[16], "B-tag of 4th jet", 13, 0., 1.)

h_QCD[list_histos[17]] = ROOT.TH1F(list_histos[17], "m_{12} invariant mass", 15, 0., 600.)
h_QCD[list_histos[18]] = ROOT.TH1F(list_histos[18], "m_{34} invariant mass", 15, 0., 600.)
h_QCD[list_histos[19]] = ROOT.TH1F(list_histos[19], "No. of primary verticies", 50, 0., 50.)

#Defining 2D Histograms/Corelation's
h_ma1_ma2 = ROOT.TH2F("h_ma1_ma2", "m_{12} vs m_{34}", 15, 0., 600., 15, 0., 600.)
h_mh_ma1  = ROOT.TH2F("h_mh_ma2", "m_{H} vs m_{34}", 20, 200., 1000., 15, 0., 600.)
h_ma1_ma2_sig = ROOT.TH2F("h_ma1_ma2_sig", "m_{12} vs m_{34}", 15, 0., 600., 15, 0., 600.)
h_mh_ma1_sig  = ROOT.TH2F("h_mh_ma2_sig", "m_{H} vs m_{34}", 20, 200., 1000., 15, 0., 600.)


##Graphics stuff
canvas = dict()
canvas_sig = dict()
for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname)
    canvas_sig[hname] = ROOT.TCanvas(hname+"_sig",hname+"_sig")
leg1 = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg1.SetHeader("Samples considered")

Nsig_passed = 0.
Nbkg_passed = 0.

##Loop on samples, and then on events, and merge QCD stuff
idx_sample = 0
for name_sample in samplename_list:

    if "Signal" in name_sample and not name_sample == myWF.sig_samplename:
        continue

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

        jet1eta = jet1_4mom.Eta()
        jet2eta = jet2_4mom.Eta()
        jet3eta = jet3_4mom.Eta()
        jet4eta = jet4_4mom.Eta()

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

        if delta_Phi > 3.14:
            delta_Phi = 6.28 - delta_Phi

        pt_pair1 = p_pair1.Pt()
        pt_pair2 = p_pair2.Pt()

        totaljets_4mom = p_pair1 + p_pair2

        m_pair1 = p_pair1.M()
        m_pair2 = p_pair2.M()
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
            if select_all_but_one("h_jet1eta"):
                h_QCD["h_jet1eta"].Fill(jet1eta,norm_factor)
            if select_all_but_one("h_jet2eta"):
                h_QCD["h_jet2eta"].Fill(jet2eta,norm_factor)
            if select_all_but_one("h_jet3eta"):
                h_QCD["h_jet3eta"].Fill(jet3eta,norm_factor)
            if select_all_but_one("h_jet4eta"):
                h_QCD["h_jet4eta"].Fill(jet4eta,norm_factor)
            if select_all_but_one("h_jet1Btag"):
                h_QCD["h_jet1Btag"].Fill(mytree.jet1Btag,norm_factor)
            if select_all_but_one("h_jet2Btag"):
                h_QCD["h_jet2Btag"].Fill(mytree.jet2Btag,norm_factor)
            if select_all_but_one("h_jet3Btag"):
                h_QCD["h_jet3Btag"].Fill(mytree.jet3Btag,norm_factor)
            if select_all_but_one("h_jet4Btag"):
                h_QCD["h_jet4Btag"].Fill(mytree.jet4Btag,norm_factor)	
            if select_all_but_one("h_m12"):
                h_QCD["h_m12"].Fill(m_pair1,norm_factor)
            if select_all_but_one("h_m34"):
                h_QCD["h_m34"].Fill(m_pair2,norm_factor)
	    if select_all_but_one("h_nPv"):
                h_QCD["h_nPv"].Fill(mytree.N_nPv,norm_factor) #added
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
            if select_all_but_one("h_jet1eta"):
                h_base[name_sample+"h_jet1eta"].Fill(jet1eta,norm_factor)
            if select_all_but_one("h_jet2eta"):
                h_base[name_sample+"h_jet2eta"].Fill(jet2eta,norm_factor)
            if select_all_but_one("h_jet3eta"):
                h_base[name_sample+"h_jet3eta"].Fill(jet3eta,norm_factor)
            if select_all_but_one("h_jet4eta"):
                h_base[name_sample+"h_jet4eta"].Fill(jet4eta,norm_factor)
            if select_all_but_one("h_jet1Btag"):
                h_base[name_sample+"h_jet1Btag"].Fill(mytree.jet1Btag,norm_factor)
            if select_all_but_one("h_jet2Btag"):
                h_base[name_sample+"h_jet2Btag"].Fill(mytree.jet2Btag,norm_factor)
            if select_all_but_one("h_jet3Btag"):
                h_base[name_sample+"h_jet3Btag"].Fill(mytree.jet3Btag,norm_factor)
            if select_all_but_one("h_jet4Btag"):
                h_base[name_sample+"h_jet4Btag"].Fill(mytree.jet4Btag,norm_factor)
            if select_all_but_one("h_m12"):
                h_base[name_sample+"h_m12"].Fill(m_pair1,norm_factor)
            if select_all_but_one("h_m34"):
                h_base[name_sample+"h_m34"].Fill(m_pair2,norm_factor)
	    if select_all_but_one("h_nPv"):
                h_base[name_sample+"h_nPv"].Fill(mytree.N_nPv,norm_factor)

        ##Now the 2D plots
        if select_all_but_one("h_mh_ma1"):
            if name_sample == myWF.sig_samplename:
                h_ma1_ma2_sig.Fill(m_pair1,m_pair2,norm_factor)
            else:
                h_ma1_ma2.Fill(m_pair1,m_pair2,norm_factor)
        if select_all_but_one("h_mh_ma1"):
            if name_sample == myWF.sig_samplename:
                h_mh_ma1_sig.Fill(m4b,m_pair1,norm_factor)   
            else:
                h_mh_ma1.Fill(m4b,m_pair1,norm_factor)   

        #Count the events
        if select_all_but_one("actually all"):
            if name_sample == myWF.sig_samplename:
                Nsig_passed += norm_factor
            else:
                Nbkg_passed += norm_factor

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
    #canvas[hname].SetLogy(1)
    hs[hname].Draw("histo")
    if signal_magnify != 1:
        h_base[myWF.sig_samplename+hname].Scale(signal_magnify)
    leg1.Draw()

    canvas[hname].SaveAs("plots/tree_" + hname + ".gif")
    #canvas[hname].SaveAs("plots/tree_" + hname + ".C")
    canvas[hname].SaveAs("plots/tree_" + hname + ".pdf")

##Signal only plots to be produced
for hname in list_histos:
    if hname == "h_m4b" or hname == "h_m12" or  hname == "h_m34":
        canvas_sig[hname].cd()

        ##Remove some style fancyness
        h_base[myWF.sig_samplename+hname].SetFillColor(0)
        h_base[myWF.sig_samplename+hname].SetLineStyle(0)

        ROOT.gStyle.SetOptStat(1)
        h_base[myWF.sig_samplename+hname].Draw("histo")
        canvas_sig[hname].SaveAs("plots/tree_" + hname + "_signal.gif")
        #canvas_sig[hname].SaveAs("plots/tree_" + hname + "_signal.C")
        canvas_sig[hname].SaveAs("plots/tree_" + hname + "_signal.pdf")

##Now the 2D plots
ROOT.gStyle.SetOptStat(0)

canvas_ma1_ma2 = ROOT.TCanvas("h_ma1_ma2","h_ma1_ma2")
canvas_ma1_ma2.cd()
h_ma1_ma2.Draw("COLZ")
canvas_ma1_ma2.SaveAs("plots/tree_h_ma1_ma2.gif")
#canvas_ma1_ma2.SaveAs("plots/tree_h_ma1_ma2.C")
canvas_ma1_ma2.SaveAs("plots/tree_h_ma1_ma2.pdf")

canvas_ma1_ma2_sig = ROOT.TCanvas("h_ma1_ma2_sig","h_ma1_ma2_sig")
canvas_ma1_ma2_sig.cd()
h_ma1_ma2_sig.Draw("COLZ")
canvas_ma1_ma2_sig.SaveAs("plots/tree_h_ma1_ma2_signal.gif")
#canvas_ma1_ma2_sig.SaveAs("plots/tree_h_ma1_ma2_signal.C")
canvas_ma1_ma2_sig.SaveAs("plots/tree_h_ma1_ma2_signal.pfd")

canvas_mh_ma1 = ROOT.TCanvas("h_mh_ma1","h_mh_ma1")
canvas_mh_ma1.cd()
h_mh_ma1.Draw("COLZ")
canvas_mh_ma1.SaveAs("plots/tree_h_mh_ma1.gif")
#canvas_mh_ma1.SaveAs("plots/tree_h_mh_ma1.C")
canvas_mh_ma1.SaveAs("plots/tree_h_mh_ma1.pdf")

canvas_mh_ma1_sig = ROOT.TCanvas("h_mh_ma1_sig","h_mh_ma1_sig")
canvas_mh_ma1_sig.cd()
h_mh_ma1_sig.Draw("COLZ")
canvas_mh_ma1_sig.SaveAs("plots/tree_h_mh_ma1_signal.gif")
#canvas_mh_ma1_sig.SaveAs("plots/tree_h_mh_ma1_signal.C")
canvas_mh_ma1_sig.SaveAs("plots/tree_h_mh_ma1_signal.pdf")

print "Number of expected events for ", luminosity_norm, " fb-1"
print "Number of signal events = ", Nsig_passed
print "Number of background events = ", Nbkg_passed
print "Significance S/sqrt(B) = ", Nsig_passed/math.sqrt(Nbkg_passed)
print "Significance S/sqrt(B + deltaB^2) = ", Nsig_passed/(math.sqrt(Nbkg_passed) + 0.2*Nbkg_passed)
print "Significance S/sqrt(S+B) = ", Nsig_passed/math.sqrt(Nsig_passed + Nbkg_passed)
print "\nAll the intresting plots have been produced..!"
message =raw_input('\nHit the Enter Key to exit the program! ')
print(message)
print "Thank You and Good Bye.!\n\n "
exit()
