import ROOT
import os
import math
import numpy as np
import sys

sys.path.append("../")

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H800_A300")

#Normalize to this luminsity, in fb-1
luminosity_norm = 36.46

output_dir = "inputs/"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#Here starts the program
Norm_Map = myWF.get_normalizations_map()

#Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

Nsig_passed = 0.
Nbkg_passed = 0.

fOut = ROOT.TFile(output_dir+"Nominal_training.root","RECREATE")

tree_Signal = ROOT.TTree("tree_signal","The signal tree")
tree_Background = ROOT.TTree("tree_background","The background tree")

j1_btag = np.zeros(1, dtype=float)
j2_btag = np.zeros(1, dtype=float)
j3_btag = np.zeros(1, dtype=float)
j4_btag = np.zeros(1, dtype=float)

j1_pt = np.zeros(1, dtype=float)
j2_pt = np.zeros(1, dtype=float)
j3_pt = np.zeros(1, dtype=float)
j4_pt = np.zeros(1, dtype=float)

delta_phi = np.zeros(1, dtype=float)
delta_eta = np.zeros(1, dtype=float)
abs_massRatio_jetpair = np.zeros(1, dtype=float)

evt_weight = np.zeros(1, dtype=float)

tree_Signal.Branch('j1_btag', j1_btag, 'j1_btag/D')
tree_Signal.Branch('j2_btag', j2_btag, 'j2_btag/D')
tree_Signal.Branch('j3_btag', j3_btag, 'j3_btag/D')
tree_Signal.Branch('j4_btag', j4_btag, 'j4_btag/D')

tree_Signal.Branch('j1_pt', j1_pt, 'j1_pt/D')
tree_Signal.Branch('j2_pt', j2_pt, 'j2_pt/D')
tree_Signal.Branch('j3_pt', j3_pt, 'j3_pt/D')
tree_Signal.Branch('j4_pt', j4_pt, 'j4_pt/D')

tree_Signal.Branch('delta_phi', delta_phi, 'delta_phi/D')
tree_Signal.Branch('delta_eta', delta_eta, 'delta_eta/D')
tree_Signal.Branch('abs_massRatio_jetpair', abs_massRatio_jetpair, 'abs_massRatio_jetpair/D')

tree_Signal.Branch('evt_weight', evt_weight, 'evt_weight/D')

tree_Background.Branch('j1_btag', j1_btag, 'j1_btag/D')
tree_Background.Branch('j2_btag', j2_btag, 'j2_btag/D')
tree_Background.Branch('j3_btag', j3_btag, 'j3_btag/D')
tree_Background.Branch('j4_btag', j4_btag, 'j4_btag/D')

tree_Background.Branch('j1_pt', j1_pt, 'j1_pt/D')
tree_Background.Branch('j2_pt', j2_pt, 'j2_pt/D')
tree_Background.Branch('j3_pt', j3_pt, 'j3_pt/D')
tree_Background.Branch('j4_pt', j4_pt, 'j4_pt/D')

tree_Background.Branch('delta_phi', delta_phi, 'delta_phi/D')
tree_Background.Branch('delta_eta', delta_eta, 'delta_eta/D')
tree_Background.Branch('abs_massRatio_jetpair', abs_massRatio_jetpair, 'abs_massRatio_jetpair/D')

tree_Background.Branch('evt_weight', evt_weight, 'evt_weight/D')

for sample_name in samplename_list:

    isSignal = False

    if "Signal" in sample_name:
        isSignal = True
        if not sample_name == myWF.sig_samplename:
            continue

    norm_factor = Norm_Map[sample_name]*luminosity_norm

    mytree = root_file[sample_name].Get("HAA4bAnalysis/mytree")               # for background only
 
    print "Processing Sample ", sample_name  #///

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break
        nb = mytree.GetEntry(jentry )
        if nb <= 0:
            continue

        PU_Weight = mytree.PUWeight
        Event_Weight = norm_factor * PU_Weight

        evt_weight[0] = Event_Weight

        jet1_4mom = mytree.jet1_4mom
        jet2_4mom = mytree.jet2_4mom
        jet3_4mom = mytree.jet3_4mom
        jet4_4mom = mytree.jet4_4mom

        jet1_4mom_fit = mytree.jet1_4mom_fit
        jet2_4mom_fit = mytree.jet2_4mom_fit
        jet3_4mom_fit = mytree.jet3_4mom_fit
        jet4_4mom_fit = mytree.jet4_4mom_fit

        j1_pt[0] = jet1_4mom.Pt()
        j2_pt[0] = jet2_4mom.Pt()
        j3_pt[0] = jet3_4mom.Pt()
        j4_pt[0] = jet4_4mom.Pt()

        if j4_pt[0] < 50.:
            continue

        j1_btag[0] = mytree.jet1Btag
        j2_btag[0] = mytree.jet2Btag
        j3_btag[0] = mytree.jet3Btag
        j4_btag[0] = mytree.jet4Btag

        if j1_btag[0] < 0.:
            print "jet1 btag = ", j1_btag[0]
            print " jet1 pt = ", j1_pt[0]

        #get best combinations of jets before fit
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

        #get best combinations of jets after fit
        combination_flag_fit = myWF.get_best_combination(jet1_4mom_fit,jet2_4mom_fit,jet3_4mom_fit,jet4_4mom_fit);
        if combination_flag_fit == 1:
            p_pair1_fit = jet1_4mom_fit + jet2_4mom_fit
            p_pair2_fit = jet3_4mom_fit + jet4_4mom_fit
        elif combination_flag_fit == 2:
            p_pair1_fit = jet1_4mom_fit + jet3_4mom_fit
            p_pair2_fit = jet2_4mom_fit + jet4_4mom_fit
        elif combination_flag_fit == 3:
            p_pair1_fit = jet1_4mom_fit + jet4_4mom_fit
            p_pair2_fit = jet2_4mom_fit + jet3_4mom_fit

        delta_Phi = abs(p_pair1.Phi() - p_pair2.Phi())
        delta_Eta = p_pair1.Eta() - p_pair2.Eta()

        if delta_Phi > 3.14:
            delta_Phi = 6.28 - delta_Phi

        delta_phi[0] = delta_Phi
        delta_eta[0] = delta_Eta

        #abs(mass_diff_pair1/mass_add_pair2) variable before fit
        diff_p_pair1 = p_pair1 - p_pair2
        add_p_pair2 = p_pair1 + p_pair2
        mass_diff_pair1 = diff_p_pair1.M()
        mass_add_pair2 = add_p_pair2.M()
        if not mass_add_pair2 == 0.:
            abs_massRatio_jetpair[0] = abs(mass_diff_pair1/mass_add_pair2)
        else:
            abs_massRatio_jetpair[0] = 0.
    
        #Count the events
        if sample_name == myWF.sig_samplename:
            Nsig_passed += norm_factor
        else:
            Nbkg_passed += norm_factor

        if isSignal:
            tree_Signal.Fill()
        else:
            tree_Background.Fill()
            

print "Finished runnning over samples!"

print "Number of expected events for ", luminosity_norm, " in fb-1"
print "Number of signal events used = ", Nsig_passed
print "Number of background events used = ", Nbkg_passed

tree_Signal.Write()
tree_Background.Write()

fOut.Close()

exit()
