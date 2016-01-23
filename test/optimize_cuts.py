import ROOT
import math
import numpy as np

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H500_A200")

def is_Event_selected(jet1_btag,jet2_btag,jet3_btag):
    if jet1_btag > 0.97 and jet2_btag > 0.97 and jet3_btag > 0.97: return True
    #if jet1_btag > 0.97 and jet2_btag > 0.97: return True
    #if jet1_btag > 0.97: return True
    return False

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

steps_cut1 = 20
steps_cut2 = 20
steps_cut3 = 20
steps_cut4 = 20

cut_Nbkg = [[[[0. for x in range(steps_cut4)] for x in range(steps_cut3)] for x in range(steps_cut2)] for x in range(steps_cut1)]
cut_Nsig = [[[[0. for x in range(steps_cut4)] for x in range(steps_cut3)] for x in range(steps_cut2)] for x in range(steps_cut1)]

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

        if not is_Event_selected(mytree.jet1Btag,mytree.jet2Btag,mytree.jet3Btag):
            continue

        jet1pt = mytree.jet1pt
        jet2pt = mytree.jet2pt
        jet3pt = mytree.jet3pt
        delta_Phi = abs(mytree.delta_Phi_pair)
        
        for icut1 in xrange(0,steps_cut1):
            cut1_value = 50. + 5.*icut1

            if jet1pt < cut1_value:
                continue

            for icut2 in xrange(0,steps_cut2):
                cut2_value = 0. + 0.1*icut2

                if delta_Phi < cut2_value:
                    continue

                for icut3 in xrange(0,steps_cut3):
                    cut3_value = 30. + 5.*icut3

                    if jet2pt < cut3_value:
                        continue

                    if cut3_value > cut1_value:
                        continue

                    for icut4 in xrange(0,steps_cut4):
                        cut4_value = 30. + 5.*icut4
                        
                        if jet3pt < cut4_value:
                            continue

                        if cut4_value > cut3_value:
                            continue

                        if name_sample == myWF.sig_samplename:
                            cut_Nsig[icut1][icut2][icut3][icut4] = cut_Nsig[icut1][icut2][icut3][icut4] + norm_factor
                        else:
                            cut_Nbkg[icut1][icut2][icut3][icut4] = cut_Nbkg[icut1][icut2][icut3][icut4] + norm_factor

##Calculate the significance
signif_list = [[[[0. for x in range(steps_cut4)] for x in range(steps_cut3)] for x in range(steps_cut2)] for x in range(steps_cut1)]
cut1_x_list = [0.]*steps_cut1
cut2_x_list = [0.]*steps_cut2
cut3_x_list = [0.]*steps_cut3
cut4_x_list = [0.]*steps_cut4

signif_max = -1.
for icut1 in range(0,steps_cut1):
    cut1_x_list[icut1] = 50. + 5.*icut1

    for icut2 in range(0,steps_cut2):
        cut2_x_list[icut2] = 0. + 0.1*icut2

        for icut3 in range(0,steps_cut3):
            cut3_x_list[icut3] = 30. + 5.*icut3

            for icut4 in range(0,steps_cut4):
                cut4_x_list[icut4] = 30. + 5.*icut4

                if cut_Nbkg[icut1][icut2][icut3][icut4] != 0:
                    signif_list[icut1][icut2][icut3][icut4] = cut_Nsig[icut1][icut2][icut3][icut4]/(math.sqrt(cut_Nbkg[icut1][icut2][icut3][icut4]) + 0.2*cut_Nbkg[icut1][icut2][icut3][icut4])
                else:
                    signif_list[icut1][icut2][icut3][icut4] = 0.

                if signif_list[icut1][icut2][icut3][icut4] > signif_max:
                    signif_max = signif_list[icut1][icut2][icut3][icut4]
                    cut1_max = icut1
                    cut2_max = icut2
                    cut3_max = icut3
                    cut4_max = icut4

print "The cut1 value is ", 50. + 5.*cut1_max
print "The cut2 value is ", 0. + 0.1*cut2_max
print "The cut3 value is ", 30. + 5.*cut3_max
print "The cut4 value is ", 30. + 5.*cut4_max
print "Max significance is ", signif_max

##Plot the significance as function of the cuts
cut1_x = np.array(cut1_x_list)
cut2_x = np.array(cut2_x_list)
cut3_x = np.array(cut3_x_list)
cut4_x = np.array(cut4_x_list)

cut1_y_list = [0.]*steps_cut1
cut2_y_list = [0.]*steps_cut2
cut3_y_list = [0.]*steps_cut3
cut4_y_list = [0.]*steps_cut4

for icut1 in range(0,steps_cut1):

    for icut2 in range(0,steps_cut2):

        for icut3 in range(0,steps_cut3):
            
            for icut4 in range(0,steps_cut4):

                if icut2 == cut2_max and icut3 == cut3_max and icut4 == cut4_max:
                    cut1_y_list[icut1] = signif_list[icut1][icut2][icut3][icut4]

                if icut1 == cut1_max and icut3 == cut3_max and icut4 == cut4_max:
                    cut2_y_list[icut2] = signif_list[icut1][icut2][icut3][icut4]

                if icut1 == cut1_max and icut2 == cut2_max and icut4 == cut4_max:
                    cut3_y_list[icut3] = signif_list[icut1][icut2][icut3][icut4]

                if icut1 == cut1_max and icut2 == cut2_max and icut3 == cut3_max:
                    cut4_y_list[icut4] = signif_list[icut1][icut2][icut3][icut4]


cut1_y = np.array(cut1_y_list)
cut2_y = np.array(cut2_y_list)
cut3_y = np.array(cut3_y_list)
cut4_y = np.array(cut4_y_list)

graph_cut1 = ROOT.TGraph(steps_cut1,cut1_x,cut1_y)
graph_cut2 = ROOT.TGraph(steps_cut2,cut2_x,cut2_y)
graph_cut3 = ROOT.TGraph(steps_cut3,cut3_x,cut3_y)
graph_cut4 = ROOT.TGraph(steps_cut4,cut4_x,cut4_y)

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
