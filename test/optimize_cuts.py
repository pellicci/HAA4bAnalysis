import ROOT
import math
import numpy as np

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H500_A200")

def is_Event_selected(jet_btag,jet_pt):
    """Save events according to some basic selection criteria"""
    btag_cut = jet_btag[0] > 0.97 and jet_btag[1] > 0.97 and jet_btag[2] > 0.97 and jet_btag[3] > 0.97
    #pt_cut = jet_pt[0] > 165. and jet_pt[1] > 130. and jet_pt[2] > 130. and jet_pt[3] > 110.
    pt_cut = jet_pt[3] > 50.

    return btag_cut and pt_cut

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

steps_cut1 = 20
cut1_init = 0.
cut1_stepsize = 10.

steps_cut2 = 25
cut2_init = 0.
cut2_stepsize = 0.1

steps_cut3 = 20
cut3_init = 0.
cut3_stepsize = 10.

steps_cut4 = 15
cut4_init = 100.
cut4_stepsize = 10.

steps_cut5 = 15
cut5_init = 50.
cut5_stepsize = 10.

steps_cut6 = 15
cut6_init = 50.
cut6_stepsize = 10.

cut_Nbkg = [[[[[[0. for x in range(steps_cut6)] for x in range(steps_cut5)] for x in range(steps_cut4)] for x in range(steps_cut3)] for x in range(steps_cut2)] for x in range(steps_cut1)]
cut_Nsig = [[[[[[0. for x in range(steps_cut6)] for x in range(steps_cut5)] for x in range(steps_cut4)] for x in range(steps_cut3)] for x in range(steps_cut2)] for x in range(steps_cut1)]

##collect all the root files
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##loop on the samples and on the cuts and calculate N_bkg and N_sig
for name_sample in samplename_list:

    norm_factor = Norm_Map[name_sample]
    mytree = root_file[name_sample].Get("HZZ4bAnalysis/mytree")

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break

        nb = mytree.GetEntry( jentry )
        if nb <= 0:
            continue

        jet_btag = [mytree.jet1Btag, mytree.jet2Btag, mytree.jet3Btag, mytree.jet4Btag]

        jet1_4mom = mytree.jet1_4mom
        jet2_4mom = mytree.jet2_4mom
        jet3_4mom = mytree.jet3_4mom
        jet4_4mom = mytree.jet4_4mom

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

        jet1pt = jet1_4mom.Pt()
        jet2pt = jet2_4mom.Pt()
        jet3pt = jet3_4mom.Pt()
        jet4pt = jet4_4mom.Pt()

        jet_pt = [jet1pt, jet2pt, jet3pt, jet4pt]

        if not is_Event_selected(jet_btag,jet_pt):
            continue

        delta_Phi = abs(p_pair1.Phi() - p_pair2.Phi())
        delta_Eta = abs(p_pair1.Eta() - p_pair2.Eta())
        
        for icut1 in xrange(steps_cut1):
            cut1_value = cut1_init + cut1_stepsize*icut1

            if p_pair1.Pt() < cut1_value:
                continue

            for icut2 in xrange(steps_cut2):
                cut2_value = cut2_init + cut2_stepsize*icut2

                if delta_Phi < cut2_value:
                    continue

                for icut3 in xrange(steps_cut3):
                    cut3_value = cut3_init + cut3_stepsize*icut3

                    if p_pair2.Pt() < cut3_value:
                        continue

                    for icut4 in xrange(steps_cut4):
                        cut4_value = cut4_init + cut4_stepsize*icut4
                        
                        if jet1pt < cut4_value:
                            continue

                        for icut5 in xrange(steps_cut5):
                            cut5_value = cut5_init + cut5_stepsize*icut5

                            if jet2pt < cut5_value:
                                continue

                            if cut5_value > cut4_value:
                                continue

                            for icut6 in xrange(steps_cut6):
                                cut6_value = cut6_init + cut6_stepsize*icut6

                                if jet3pt < cut6_value:
                                    continue

                                if cut6_value > cut5_value:
                                    continue
                            
                                if name_sample == myWF.sig_samplename:
                                    cut_Nsig[icut1][icut2][icut3][icut4][icut5][icut6] = cut_Nsig[icut1][icut2][icut3][icut4][icut5][icut6] + norm_factor
                                else:
                                    cut_Nbkg[icut1][icut2][icut3][icut4][icut5][icut6] = cut_Nbkg[icut1][icut2][icut3][icut4][icut5][icut6] + norm_factor

##Calculate the significance
signif_list = [[[[[[0. for x in range(steps_cut6)] for x in range(steps_cut5)] for x in range(steps_cut4)] for x in range(steps_cut3)] for x in range(steps_cut2)] for x in range(steps_cut1)]

cut1_x_list = [cut1_init + cut1_stepsize*icut for icut in range(steps_cut1)]
cut2_x_list = [cut2_init + cut2_stepsize*icut for icut in range(steps_cut2)]
cut3_x_list = [cut3_init + cut3_stepsize*icut for icut in range(steps_cut3)]
cut4_x_list = [cut4_init + cut4_stepsize*icut for icut in range(steps_cut4)]
cut5_x_list = [cut5_init + cut5_stepsize*icut for icut in range(steps_cut5)]
cut6_x_list = [cut6_init + cut6_stepsize*icut for icut in range(steps_cut6)]

signif_max = -1.
for icut1 in xrange(steps_cut1):
    for icut2 in xrange(steps_cut2):
        for icut3 in xrange(steps_cut3):
            for icut4 in xrange(steps_cut4):
                for icut5 in xrange(steps_cut5):
                    for icut6 in xrange(steps_cut6):

                        if cut_Nbkg[icut1][icut2][icut3][icut4][icut5][icut6] != 0:
                            signif_list[icut1][icut2][icut3][icut4][icut5][icut6] = cut_Nsig[icut1][icut2][icut3][icut4][icut5][icut6]/(
                                math.sqrt(cut_Nbkg[icut1][icut2][icut3][icut4][icut5][icut6]) + 0.2*cut_Nbkg[icut1][icut2][icut3][icut4][icut5][icut6])
                        else:
                            signif_list[icut1][icut2][icut3][icut4][icut5][icut6] = 0.

                        if signif_list[icut1][icut2][icut3][icut4][icut5][icut6] > signif_max:
                            signif_max = signif_list[icut1][icut2][icut3][icut4][icut5][icut6] 
                            cut1_max = icut1
                            cut2_max = icut2
                            cut3_max = icut3
                            cut4_max = icut4
                            cut5_max = icut5
                            cut6_max = icut6

print "The cut1 value is ", cut1_init + cut1_stepsize*cut1_max
print "The cut2 value is ", cut2_init + cut2_stepsize*cut2_max
print "The cut3 value is ", cut3_init + cut3_stepsize*cut3_max
print "The cut4 value is ", cut4_init + cut4_stepsize*cut4_max
print "The cut5 value is ", cut5_init + cut5_stepsize*cut5_max
print "The cut6 value is ", cut6_init + cut6_stepsize*cut6_max
print "Max significance is ", signif_max

##Plot the significance as function of the cuts
cut1_x = np.array(cut1_x_list)
cut2_x = np.array(cut2_x_list)
cut3_x = np.array(cut3_x_list)
cut4_x = np.array(cut4_x_list)
cut5_x = np.array(cut5_x_list)
cut6_x = np.array(cut6_x_list)

cut1_y_list = [signif_list[icut][cut2_max][cut3_max][cut4_max][cut5_max][cut6_max] for icut in xrange(steps_cut1)]
cut2_y_list = [signif_list[cut1_max][icut][cut3_max][cut4_max][cut5_max][cut6_max] for icut in xrange(steps_cut2)]
cut3_y_list = [signif_list[cut1_max][cut2_max][icut][cut4_max][cut5_max][cut6_max] for icut in xrange(steps_cut3)]
cut4_y_list = [signif_list[cut1_max][cut2_max][cut3_max][icut][cut5_max][cut6_max] for icut in xrange(steps_cut4)]
cut5_y_list = [signif_list[cut1_max][cut2_max][cut3_max][cut4_max][icut][cut6_max] for icut in xrange(steps_cut5)]
cut6_y_list = [signif_list[cut1_max][cut2_max][cut3_max][cut4_max][cut5_max][icut] for icut in xrange(steps_cut6)]

cut1_y = np.array(cut1_y_list)
cut2_y = np.array(cut2_y_list)
cut3_y = np.array(cut3_y_list)
cut4_y = np.array(cut4_y_list)
cut5_y = np.array(cut5_y_list)
cut6_y = np.array(cut6_y_list)

graph_cut1 = ROOT.TGraph(steps_cut1,cut1_x,cut1_y)
graph_cut2 = ROOT.TGraph(steps_cut2,cut2_x,cut2_y)
graph_cut3 = ROOT.TGraph(steps_cut3,cut3_x,cut3_y)
graph_cut4 = ROOT.TGraph(steps_cut4,cut4_x,cut4_y)
graph_cut5 = ROOT.TGraph(steps_cut5,cut5_x,cut5_y)
graph_cut6 = ROOT.TGraph(steps_cut6,cut6_x,cut6_y)

c1 = ROOT.TCanvas("c1","c1")
c1.cd()
graph_cut1.Draw("AC*")
c1.SaveAs("plots/cut1_signif.gif")

c2 = ROOT.TCanvas("c2","c2")
c2.cd()
graph_cut2.Draw("AC*")
c2.SaveAs("plots/cut2_signif.gif")

c3 = ROOT.TCanvas("c3","c3")
c3.cd()
graph_cut3.Draw("AC*")
c3.SaveAs("plots/cut3_signif.gif")

c4 = ROOT.TCanvas("c4","c4")
c4.cd()
graph_cut4.Draw("AC*")
c4.SaveAs("plots/cut4_signif.gif")

c5 = ROOT.TCanvas("c5","c5")
c5.cd()
graph_cut5.Draw("AC*")
c5.SaveAs("plots/cut5_signif.gif")

c6 = ROOT.TCanvas("c6","c6")
c6.cd()
graph_cut6.Draw("AC*")
c6.SaveAs("plots/cut6_signif.gif")
