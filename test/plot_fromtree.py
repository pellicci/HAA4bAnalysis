import ROOT
import os
import math
#import numpy as np

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H800_A300")

##Global constants

##This is for mH = 800 and mA = 300
#PT1_MIN = 150.
#PT2_MIN = 150.
#PT3_MIN = 80.
#PT4_MIN = 80.
#DELTA_PHI_MIN = 2.4
#DELTA_ETA_MAX = 10.
#PT_PAIR1_MIN = 0.
#PT_PAIR2_MIN = 0.
#TTTT

##This is for mH = 800 and mA = 300
PT1_MIN = 50.
PT2_MIN = 50.
PT3_MIN = 50.
PT4_MIN = 50.
ETA1_MAX = 5.
ETA2_MAX = 5.
ETA3_MAX = 5.
ETA4_MAX = 5.
DELTA_PHI_MIN = 0.
DELTA_ETA_MAX = 10.
PT_PAIR1_MIN = 0.
PT_PAIR2_MIN = 0.

JET1_BTAG = 0.80 
JET2_BTAG = 0.80 
JET3_BTAG = 0.80 
JET4_BTAG = 0.80 

##Normalize to this luminsity, in fb-1
luminosity_norm = 2.178

##Make signal histos larger
signal_magnify = 100.

output_dir = "plots"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

list_histos = ["h_jet1pt", "h_jet2pt", "h_jet3pt", "h_jet4pt", "h_delta_Phi_pair", "h_delta_Eta_pair", "h_pt_pair1", "h_pt_pair2", "h_m4b", "h_m4b_fitted", "h_jet1eta", "h_jet2eta", "h_jet3eta", "h_jet4eta", "h_jet1Btag", "h_jet2Btag", "h_jet3Btag", "h_jet4Btag", "h_m12", "h_m34","h_nPv" , "h_abs_massRatio_jetpair", "h_abs_massRatio_jetpair_afterfit"  ] #"h_njetE" 

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

##Get data file names and the corresponding rootfiles
dataName_list = myWF.get_dataSample_names()
root_file_data = myWF.get_data_root_files()

#Combining two input lists i,e sample monte carlo and data
combined_list = dict()
combined_list = samplename_list + dataName_list

#Store root file in a common list
rootfiles_combined_list = dict()
rootfiles_combined_list = dict(root_file, **root_file_data)

##Get the handlers for all the histos
hs = dict()
hs_data = dict()
h_sum_mc=dict()
h_sum_qcd=dict()
h_sum_data=dict()
h_sum_total=dict()      #total sum of mc and QCD
h_ratio=dict()
h_QCD = dict()
h_base = dict()
h_Signal = dict()

for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")
    hs_data[hname] = ROOT.THStack("hs_"+ hname, "")
##Define the histos to be created

#for sample_name in samplename_list:   #for background only
#for sample_name in dataName_list:     #for data only
for sample_name in combined_list:      #for data + background

    if "QCD" in sample_name:
        continue
    if "Signal" in sample_name and not sample_name == myWF.sig_samplename:
        continue

    h_base[sample_name+list_histos[0]]  = ROOT.TH1F(sample_name+list_histos[0], "p_{T} of the 1st jet", 25, PT1_MIN, 390.)
    h_base[sample_name+list_histos[1]]  = ROOT.TH1F(sample_name+list_histos[1], "p_{T} of the 2nd jet", 25, PT2_MIN, 400.)
    h_base[sample_name+list_histos[2]]  = ROOT.TH1F(sample_name+list_histos[2], "p_{T} of the 3rd jet", 25, PT3_MIN, 400.)
    h_base[sample_name+list_histos[3]]  = ROOT.TH1F(sample_name+list_histos[3], "p_{T} of the 4th jet", 25, PT4_MIN, 390.)
    h_base[sample_name+list_histos[4]]  = ROOT.TH1F(sample_name+list_histos[4], "#Delta#phi of the two jet pairs", 20, 0., 3.14)
    h_base[sample_name+list_histos[5]]  = ROOT.TH1F(sample_name+list_histos[5], "#Delta#eta of the two jet pairs", 25, -5., 5.)
    h_base[sample_name+list_histos[6]]  = ROOT.TH1F(sample_name+list_histos[6], "p_{T} of the first jet pair", 20, 0., 500.)
    h_base[sample_name+list_histos[7]]  = ROOT.TH1F(sample_name+list_histos[7], "p_{T} of the second jet pair", 20, 0., 500.)
    h_base[sample_name+list_histos[8]]  = ROOT.TH1F(sample_name+list_histos[8], "4jets invariant mass", 25, 100., 1200.)
    h_base[sample_name+list_histos[9]]  = ROOT.TH1F(sample_name+list_histos[9], "Fitted 4jets invariant mass", 25, 100., 1200.)
    h_base[sample_name+list_histos[10]] = ROOT.TH1F(sample_name+list_histos[10], "#eta of the 1st jet", 25, -3.5, 3.5)
    h_base[sample_name+list_histos[11]] = ROOT.TH1F(sample_name+list_histos[11], "#eta of the 2nd jet", 25, -3.5, 3.5)
    h_base[sample_name+list_histos[12]] = ROOT.TH1F(sample_name+list_histos[12], "#eta of the 3rd jet", 25, -3.5, 3.5)
    h_base[sample_name+list_histos[13]] = ROOT.TH1F(sample_name+list_histos[13], "#eta of the 4th jet", 25, -3.5, 3.5)
    h_base[sample_name+list_histos[14]] = ROOT.TH1F(sample_name+list_histos[14], "B-tag of the 1st jet", 25, JET1_BTAG+0.09, 1.)
    h_base[sample_name+list_histos[15]] = ROOT.TH1F(sample_name+list_histos[15], "B-tag of the 2nd jet", 25, JET2_BTAG+0.04, 1.)
    h_base[sample_name+list_histos[16]] = ROOT.TH1F(sample_name+list_histos[16], "B-tag of the 3rd jet", 25, JET3_BTAG, 1.)
    h_base[sample_name+list_histos[17]] = ROOT.TH1F(sample_name+list_histos[17], "B-tag of the 4th jet", 25, JET4_BTAG, 1.)
    h_base[sample_name+list_histos[18]] = ROOT.TH1F(sample_name+list_histos[18], "Invariant mass m_{12} of the first jet pair", 20, 0., 550.)
    h_base[sample_name+list_histos[19]] = ROOT.TH1F(sample_name+list_histos[19], "Invariant mass m_{34} of the second jet pair", 20, 0., 550.)
    h_base[sample_name+list_histos[20]] = ROOT.TH1F(sample_name+list_histos[20], "No. of primary verticies",25 , 0., 25.)
    h_base[sample_name+list_histos[21]] = ROOT.TH1F(sample_name+list_histos[21], "abs([M(b1b2)-M(b3b4)]/[M(b1b2)+M(b1b2])",15 , 0., 1.1)
    h_base[sample_name+list_histos[22]] = ROOT.TH1F(sample_name+list_histos[22], "abs([M(b1b2)-M(b3b4)]/[M(b1b2)+M(b1b2])",15 , 0., 1.1)
    #h_base[sample_name+list_histos[20]] = ROOT.TH1F(sample_name+list_histos[20], "No. of jet entries", 10, 0., 10.)

h_QCD[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{T} of the 1st jet", 25, PT1_MIN, 390.)
h_QCD[list_histos[1]]  = ROOT.TH1F(list_histos[1], "p_{T} of the 2nd jet", 25, PT2_MIN, 400.)
h_QCD[list_histos[2]]  = ROOT.TH1F(list_histos[2], "p_{T} of the 3rd jet", 25, PT3_MIN, 400.)
h_QCD[list_histos[3]]  = ROOT.TH1F(list_histos[3], "p_{T} of the 4th jet", 25, PT4_MIN, 390.)
h_QCD[list_histos[4]]  = ROOT.TH1F(list_histos[4], "#Delta#phi of the two jet pairs", 20, 0., 3.14)
h_QCD[list_histos[5]]  = ROOT.TH1F(list_histos[5], "#Delta#eta of the two jet pairs", 25, -5., 5.)
h_QCD[list_histos[6]]  = ROOT.TH1F(list_histos[6], "p_{T} of the first jet pair", 20, 0., 500.)
h_QCD[list_histos[7]]  = ROOT.TH1F(list_histos[7], "p_{T} of the second jet pair", 20, 0., 500.)
h_QCD[list_histos[8]]  = ROOT.TH1F(list_histos[8], "4jets invariant mass", 25, 100., 1200.)
h_QCD[list_histos[9]]  = ROOT.TH1F(list_histos[9], "Fitted 4jets invariant mass", 25, 100., 1200.)
h_QCD[list_histos[10]] = ROOT.TH1F(list_histos[10], "#eta of 1st jet", 25, -3.5, 3.5)
h_QCD[list_histos[11]] = ROOT.TH1F(list_histos[11], "#eta of 2nd jet", 25, -3.5, 3.5)
h_QCD[list_histos[12]] = ROOT.TH1F(list_histos[12], "#eta of 3rd jet", 25, -3.5, 3.5)
h_QCD[list_histos[13]] = ROOT.TH1F(list_histos[13], "#eta of 4th jet", 25, -3.5, 3.5)
h_QCD[list_histos[14]] = ROOT.TH1F(list_histos[14], "B-tag of 1st jet", 25, JET1_BTAG+0.09, 1.)
h_QCD[list_histos[15]] = ROOT.TH1F(list_histos[15], "B-tag of 2nd jet", 25, JET2_BTAG+0.04, 1.)
h_QCD[list_histos[16]] = ROOT.TH1F(list_histos[16], "B-tag of 3rd jet", 25, JET3_BTAG, 1.)
h_QCD[list_histos[17]] = ROOT.TH1F(list_histos[17], "B-tag of 4th jet", 25, JET4_BTAG, 1.)
h_QCD[list_histos[18]] = ROOT.TH1F(list_histos[18], "Invariant mass m_{12} of the first jet pair", 20, 0., 550.)
h_QCD[list_histos[19]] = ROOT.TH1F(list_histos[19], "Invariant mass m_{34} of the second jet pair", 20, 0., 550.)
h_QCD[list_histos[20]] = ROOT.TH1F(list_histos[20], "No. of primary verticies", 25, 0., 25.)
h_QCD[list_histos[21]] = ROOT.TH1F(list_histos[21], "abs([M(b1b2)-M(b3b4)]/[M(b1b2)+M(b1b2])",15 , 0., 1.1)
h_QCD[list_histos[22]] = ROOT.TH1F(list_histos[22], "abs([M(b1b2)-M(b3b4)]/[M(b1b2)+M(b1b2])",15 , 0., 1.1)


#Defining 2D Histograms/Correlation's
h_ma1_ma2 = ROOT.TH2F("h_ma1_ma2", "m_{12} vs m_{34}", 15, 0., 600., 15, 0., 600.)
h_mh_ma1  = ROOT.TH2F("h_mh_ma2", "m_{H} vs m_{34}", 20, 200., 1000., 15, 0., 600.)
h_ma1_ma2_sig = ROOT.TH2F("h_ma1_ma2_sig", "m_{12} vs m_{34}", 15, 0., 600., 15, 0., 600.)
h_mh_ma1_sig  = ROOT.TH2F("h_mh_ma2_sig", "m_{H} vs m_{34}", 20, 200., 1000., 15, 0., 600.)

##Graphics stuff
canvas = dict()
canvas_sig = dict()
pad1 = dict()
pad2 = dict()

for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname,200,50,600, 600)
    canvas_sig[hname] = ROOT.TCanvas(hname+"_sig",hname+"_sig", 200,50, 600, 600)

    pad1[hname] = ROOT.TPad(hname, hname, 0,0.288,0.9984,0.9972)   #changes from here
    pad2[hname] = ROOT.TPad(hname, hname, 0,0.004,0.9983,0.3321)    
    
leg1 = ROOT.TLegend(0.6,0.55,0.88,0.88)
leg1.SetHeader("Data and MC Samples")
leg1.SetFillColor(0)
leg1.SetBorderSize(0)

Nsig_passed = 0.
Nbkg_passed = 0.
Ndata_passed = 0.

##Loop on samples, and then on events, and merge QCD stuff
idx_sample = 0

#for name_sample in samplename_list: # for only background
#for  name_sample in dataName_list:  # for only data
for  name_sample in combined_list:    # for data + background

    if "Signal" in name_sample and not name_sample == myWF.sig_samplename:
        continue

    if not "QCD" in name_sample:
        QCDflag = False
    else:
        QCDflag = True
     
    norm_factor = Norm_Map[name_sample]*luminosity_norm

    #mytree = root_file[name_sample].Get("HAA4bAnalysis/mytree")               # for background only
    #mytree = root_file_data[name_sample].Get("HAA4bAnalysis/mytree")          #for data only
    mytree = rootfiles_combined_list[name_sample].Get("HAA4bAnalysis/mytree")  #for data + background
 
    print "Processing Sample ", name_sample  #///
    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break
        nb = mytree.GetEntry(jentry )
        if nb <= 0:
            continue
        if name_sample == "BTagCSV_D" or name_sample == "BTagCSV_C":
            norm_factor_data = 1.0
            norm_factor = norm_factor_data                                        #to take care for data normalization
            PU_Weight = 1.0 
            Event_Weight = norm_factor * PU_Weight                               #Event_Weight = 1.0 for data
            
        elif name_sample == "Signal_H800_A300":
             PU_Weight = 1.0            #temporary, because there is no pileUp information in signal sample right now.
             Event_Weight = norm_factor * PU_Weight      
             #continue  
        else:
            PU_Weight = mytree.PUWeight
            Event_Weight = norm_factor * PU_Weight
        jet1_4mom = mytree.jet1_4mom
        jet2_4mom = mytree.jet2_4mom
        jet3_4mom = mytree.jet3_4mom
        jet4_4mom = mytree.jet4_4mom

        jet1_4mom_fit = mytree.jet1_4mom_fit
        jet2_4mom_fit = mytree.jet2_4mom_fit
        jet3_4mom_fit = mytree.jet3_4mom_fit
        jet4_4mom_fit = mytree.jet4_4mom_fit

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

            p_pair1_fit = jet1_4mom_fit + jet2_4mom_fit
            p_pair2_fit = jet3_4mom_fit + jet4_4mom_fit

        elif combination_flag == 2:
            p_pair1 = jet1_4mom + jet3_4mom
            p_pair2 = jet2_4mom + jet4_4mom

            p_pair1_fit = jet1_4mom_fit + jet3_4mom_fit
            p_pair2_fit = jet2_4mom_fit + jet4_4mom_fit

        elif combination_flag == 3:
            p_pair1 = jet1_4mom + jet4_4mom
            p_pair2 = jet2_4mom + jet3_4mom

            p_pair1_fit = jet1_4mom_fit + jet4_4mom_fit
            p_pair2_fit = jet2_4mom_fit + jet3_4mom_fit

        delta_Phi = abs(p_pair1.Phi() - p_pair2.Phi())
        delta_Eta = p_pair1.Eta() - p_pair2.Eta()

        if delta_Phi > 3.14:
            delta_Phi = 6.28 - delta_Phi

        pt_pair1 = p_pair1.Pt()
        pt_pair2 = p_pair2.Pt()

        totaljets_4mom = p_pair1 + p_pair2
        m4b = totaljets_4mom.M()

        pt_pair1_fit = p_pair1_fit.Pt()
        pt_pair2_fit = p_pair2_fit.Pt()


        #pt_pair1 = p_pair1.Pt()
        #pt_pair2 = p_pair2.Pt()

        totaljets_4mom_fitted = p_pair1_fit + p_pair2_fit
        m4b_fitted = totaljets_4mom_fitted.M()

        #abs(mass_diff_pair1/mass_add_pair2) variable before fit
        diff_p_pair1 = p_pair1 - p_pair2
        add_p_pair2 = p_pair1 + p_pair2
        mass_diff_pair1 = diff_p_pair1.M()
        mass_add_pair2 = add_p_pair2.M()

        if not mass_add_pair2==0:
           abs_massRatio_jetpair = abs(mass_diff_pair1/mass_add_pair2)
    
        #abs(mass_diff_pair1/mass_add_pair2) variable after fit
        diff_p_pair1_fit = p_pair1_fit - p_pair2_fit
        add_p_pair2_fit = p_pair1_fit + p_pair2_fit
        mass_diff_pair1_fit = diff_p_pair1_fit.M()
        mass_add_pair2_fit = add_p_pair2_fit.M()

        if not mass_add_pair2_fit==0:
           abs_massRatio_jetpair_afterfit = abs(mass_diff_pair1_fit/mass_add_pair2_fit)

        #p_pair1_fit = jet1_4mom_fit + jet3_4mom_fit
        #p_pair2_fit = jet2_4mom_fit + jet4_4mom_fit



        if QCDflag:
            if select_all_but_one("h_jet1pt"):
                h_QCD["h_jet1pt"].Fill(jet1pt,Event_Weight)
            if select_all_but_one("h_jet2pt"):
                h_QCD["h_jet2pt"].Fill(jet2pt,Event_Weight)
            if select_all_but_one("h_jet3pt"):
                h_QCD["h_jet3pt"].Fill(jet3pt,Event_Weight)
            if select_all_but_one("h_jet4pt"):
                h_QCD["h_jet4pt"].Fill(jet4pt,Event_Weight)
            if select_all_but_one("h_delta_Phi_pair"):
                h_QCD["h_delta_Phi_pair"].Fill(delta_Phi,Event_Weight)
            if select_all_but_one("h_delta_Eta_pair"):
                h_QCD["h_delta_Eta_pair"].Fill(delta_Eta,Event_Weight)
            if select_all_but_one("h_pt_pair1"):
                h_QCD["h_pt_pair1"].Fill(pt_pair1,Event_Weight)
            if select_all_but_one("h_pt_pair2"):
                h_QCD["h_pt_pair2"].Fill(pt_pair2,Event_Weight)
            if select_all_but_one("h_m4b"):
                h_QCD["h_m4b"].Fill(m4b,Event_Weight)
            if select_all_but_one("h_m4b_fitted"):
                h_QCD["h_m4b_fitted"].Fill(m4b_fitted,Event_Weight)
            if select_all_but_one("h_jet1eta"):
                h_QCD["h_jet1eta"].Fill(jet1eta,Event_Weight)
            if select_all_but_one("h_jet2eta"):
                h_QCD["h_jet2eta"].Fill(jet2eta,Event_Weight)
            if select_all_but_one("h_jet3eta"):
                h_QCD["h_jet3eta"].Fill(jet3eta,Event_Weight)
            if select_all_but_one("h_jet4eta"):
                h_QCD["h_jet4eta"].Fill(jet4eta,Event_Weight)
            if select_all_but_one("h_jet1Btag"):
                h_QCD["h_jet1Btag"].Fill(mytree.jet1Btag,Event_Weight)
            if select_all_but_one("h_jet2Btag"):
                h_QCD["h_jet2Btag"].Fill(mytree.jet2Btag,Event_Weight)
            if select_all_but_one("h_jet3Btag"):
                h_QCD["h_jet3Btag"].Fill(mytree.jet3Btag,Event_Weight)
            if select_all_but_one("h_jet4Btag"):
                h_QCD["h_jet4Btag"].Fill(mytree.jet4Btag,Event_Weight)	
            if select_all_but_one("h_m12"):
                h_QCD["h_m12"].Fill(p_pair1.M(),Event_Weight)
            if select_all_but_one("h_m34"):
                h_QCD["h_m34"].Fill(p_pair2.M(),Event_Weight)
	    if select_all_but_one("h_nPv"):
                h_QCD["h_nPv"].Fill(mytree.N_nPv,Event_Weight)
            if select_all_but_one("h_abs_massRatio_jetpair"):
                h_QCD["h_abs_massRatio_jetpair"].Fill(abs_massRatio_jetpair,Event_Weight)
            if select_all_but_one("h_abs_massRatio_jetpair_afterfit"):
                h_QCD["h_abs_massRatio_jetpair_afterfit"].Fill(abs_massRatio_jetpair_afterfit,Event_Weight)   

        else:
            if select_all_but_one("h_jet1pt"):
                h_base[name_sample+"h_jet1pt"].Fill(jet1pt,Event_Weight)
            if select_all_but_one("h_jet2pt"):
                h_base[name_sample+"h_jet2pt"].Fill(jet2pt,Event_Weight)
            if select_all_but_one("h_jet3pt"):
                h_base[name_sample+"h_jet3pt"].Fill(jet3pt,Event_Weight)
            if select_all_but_one("h_jet4pt"):
                h_base[name_sample+"h_jet4pt"].Fill(jet4pt,Event_Weight)
            if select_all_but_one("h_delta_Phi_pair"):
                h_base[name_sample+"h_delta_Phi_pair"].Fill(delta_Phi,Event_Weight)
            if select_all_but_one("h_delta_Eta_pair"):
                h_base[name_sample+"h_delta_Eta_pair"].Fill(delta_Eta,Event_Weight)
            if select_all_but_one("h_pt_pair1"):
                h_base[name_sample+"h_pt_pair1"].Fill(pt_pair1,Event_Weight)
            if select_all_but_one("h_pt_pair2"):
                h_base[name_sample+"h_pt_pair2"].Fill(pt_pair2,Event_Weight)
            if select_all_but_one("h_m4b"):
                h_base[name_sample+"h_m4b"].Fill(m4b,Event_Weight)
            if select_all_but_one("h_m4b_fitted"):
                h_base[name_sample+"h_m4b_fitted"].Fill(m4b_fitted,Event_Weight)
            if select_all_but_one("h_jet1eta"):
                h_base[name_sample+"h_jet1eta"].Fill(jet1eta,Event_Weight)
            if select_all_but_one("h_jet2eta"):
                h_base[name_sample+"h_jet2eta"].Fill(jet2eta,Event_Weight)
            if select_all_but_one("h_jet3eta"):
                h_base[name_sample+"h_jet3eta"].Fill(jet3eta,Event_Weight)
            if select_all_but_one("h_jet4eta"):
                h_base[name_sample+"h_jet4eta"].Fill(jet4eta,Event_Weight)
            if select_all_but_one("h_jet1Btag"):
                h_base[name_sample+"h_jet1Btag"].Fill(mytree.jet1Btag,Event_Weight)
            if select_all_but_one("h_jet2Btag"):
                h_base[name_sample+"h_jet2Btag"].Fill(mytree.jet2Btag,Event_Weight)
            if select_all_but_one("h_jet3Btag"):
                h_base[name_sample+"h_jet3Btag"].Fill(mytree.jet3Btag,Event_Weight)
            if select_all_but_one("h_jet4Btag"):
                h_base[name_sample+"h_jet4Btag"].Fill(mytree.jet4Btag,Event_Weight)
            if select_all_but_one("h_m12"): 
                h_base[name_sample+"h_m12"].Fill(p_pair1.M(),Event_Weight)
            if select_all_but_one("h_m34"):
                h_base[name_sample+"h_m34"].Fill(p_pair2.M(),Event_Weight)
	    if select_all_but_one("h_nPv"):
                h_base[name_sample+"h_nPv"].Fill(mytree.N_nPv,Event_Weight)
	    if select_all_but_one("h_abs_massRatio_jetpair"):
                h_base[name_sample+"h_abs_massRatio_jetpair"].Fill(abs_massRatio_jetpair,Event_Weight)
	    if select_all_but_one("h_abs_massRatio_jetpair_afterfit"):
                h_base[name_sample+"h_abs_massRatio_jetpair_afterfit"].Fill(abs_massRatio_jetpair_afterfit, Event_Weight)


        ##Now the 2D plots
        if select_all_but_one("h_mh_ma1"):
            if name_sample == myWF.sig_samplename:
                h_ma1_ma2_sig.Fill(p_pair1.M(),p_pair2.M(),Event_Weight)
            else:
                h_ma1_ma2.Fill(p_pair1.M(),p_pair2.M(),Event_Weight)
        if select_all_but_one("h_mh_ma1"):
            if name_sample == myWF.sig_samplename:
                h_mh_ma1_sig.Fill(m4b,p_pair1.M(),Event_Weight)   
            else:
                h_mh_ma1.Fill(m4b,p_pair1.M(),Event_Weight)   

        #Count the events
        if select_all_but_one("actually all"):
            if name_sample == myWF.sig_samplename:
                Nsig_passed += norm_factor
            elif name_sample == "BTagCSV_D" or name_sample == "BTagCSV_C":
                Ndata_passed += norm_factor
            else:
                Nbkg_passed += norm_factor

    if QCDflag:
        continue

    for idx_histo,hname in enumerate(list_histos):
        h_base[name_sample+hname].SetFillColor(idx_sample+3)
        str1 = "BTagCSV_D"
        str2 = "BTagCSV_C"
        str3 = "Data"
        if name_sample ==  str1 or name_sample == str2:
                 h_base[name_sample+hname].SetMarkerStyle(20)
                 h_base[name_sample+hname].SetMarkerColor(1)
                 h_base[name_sample+hname].SetMarkerSize(1)
        if name_sample == myWF.sig_samplename:
                h_base[name_sample+hname].SetLineStyle(2)   #dashed
                h_base[name_sample+hname].SetLineColor(4)   #blue
                h_base[name_sample+hname].SetLineWidth(4)   #kind of thik

        if idx_histo == 0:
            if name_sample == str1 or name_sample ==str2: 
               if name_sample ==str1 :
                  leg1.AddEntry(h_base[name_sample+hname],str3,"p")   
            elif name_sample == myWF.sig_samplename:
                 str4 = "100 x " + name_sample
                 leg1.AddEntry(h_base[name_sample+hname],str4,"f")  #To comment when signal is has to be excluded.
            else:
                 leg1.AddEntry(h_base[name_sample+hname],name_sample,"f")

        if name_sample == str1 or name_sample == str2:
           hs_data[hname].Add(h_base[name_sample+hname])
           h_sum_data[hname]=ROOT.TH1F(h_base[name_sample+hname])
           h_sum_data[hname].Add(h_base[name_sample+hname])
                     
        #elif name_sample == myWF.sig_samplename:    #To uncomment it when the signal needs to be excluded
        #     continue

        else:
           hs[hname].Add(h_base[name_sample+hname])
           h_sum_mc[hname]=ROOT.TH1F(h_base[name_sample+hname])
           h_sum_mc[hname].Add(h_base[name_sample+hname])        

    idx_sample += 1
    if idx_sample == 1:
        idx_sample += 1

for idx_histo,hname in enumerate(list_histos):
    h_QCD[hname].SetFillColor(2) #Red
    if idx_histo == 0:
        leg1.AddEntry(h_QCD[hname],"QCD","f")
    hs[hname].Add(h_QCD[hname])
    h_sum_qcd[hname]=ROOT.TH1F(h_QCD[hname]) 
    h_sum_qcd[hname].Add(h_QCD[hname])

for hname in list_histos: 
    canvas[hname].cd()
    #canvas[hname].SetLogy(1)

    #Draw two pads pad1 and pad2
    pad1[hname].Draw()                      
    pad2[hname].Draw()
    pad1[hname].SetFillColor(0)
    pad2[hname].SetFillColor(0)
    pad1[hname].SetFrameFillColor(0)
    pad2[hname].SetFrameFillColor(0)
    pad2[hname].SetTopMargin(0.02069)
    pad1[hname].SetBottomMargin(0.10)  
    pad2[hname].SetBottomMargin(0.26)
    pad1[hname].SetFrameBorderMode(0)
    pad1[hname].SetFrameBorderMode(0)

    #Switching to first pad in the same canvas
    pad1[hname].cd()                          
    hs[hname].Draw("histo")

    #Draw the data separately
    hs_data[hname].Draw("same e1p") 
    hmax = hs[hname].GetMaximum()
    hmax_data = hs_data[hname].GetMaximum()
    if hmax>hmax_data:
        hmax= hmax
    else:
        hmax=hmax_data
    hs[hname].SetMaximum(200 + hmax)
    hs_data[hname].SetMaximum(200 + hmax)
      
    #Grapic names
    hs[hname].GetXaxis().SetTickLength(0.03) 
    hs[hname].GetXaxis().SetLabelOffset(0.006)
    hs[hname].GetYaxis().SetLabelOffset(0.007)
    hs[hname].GetYaxis().SetLabelSize(0.048)
    hs_data[hname].GetYaxis().SetLabelSize(0.048)
    hs[hname].GetYaxis().SetNdivisions(506)
    hs[hname].GetYaxis().SetTitle("Events")
    hs[hname].GetYaxis().SetTitleSize(0.045)
    hs[hname].GetXaxis().SetLabelFont(0)
    hs[hname].GetYaxis().CenterTitle()
    hs_data[hname].GetYaxis().SetNdivisions(506)

    #Switching to second pad on the same canvas
    pad2[hname].cd() 
    ROOT.gStyle.SetOptStat(0)
    h_sum_total[hname]=ROOT.TH1F(h_sum_mc[hname])
    h_sum_total[hname].Add(h_sum_qcd[hname])

    #Take the ratio between data and MC and draw it
    h_ratio[hname]= ROOT.TH1F(h_sum_data[hname])
    h_ratio[hname].Divide(h_sum_total[hname])                      
    h_ratio[hname].Draw("elp")

    #Graphic names 
    h_ratio[hname].GetXaxis().SetTitle(h_base[name_sample+hname].GetTitle())
    h_ratio[hname].SetTitle("")
    h_ratio[hname].GetXaxis().SetTickLength(0.04)    
    h_ratio[hname].GetXaxis().SetLabelFont(0)     
    h_ratio[hname].GetYaxis().SetNdivisions(506)     
    h_ratio[hname].GetYaxis().SetLabelSize(0.088)   
    h_ratio[hname].GetXaxis().SetLabelSize(0.088)     
    h_ratio[hname].GetXaxis().SetTickLength(0.09)
    h_ratio[hname].GetYaxis().SetTitle("Data/MC")
    h_ratio[hname].GetYaxis().SetTitleSize(0.08)
    h_ratio[hname].GetYaxis().SetTitleOffset(0.55)
    h_ratio[hname].GetXaxis().SetTitleOffset(1.077)
    h_ratio[hname].GetXaxis().SetTitleSize(0.097)
    ROOT.gPad.SetGridy(1)
    h_ratio[hname].GetXaxis().CenterTitle()
    h_ratio[hname].GetYaxis().CenterTitle()
    h_ratio[hname].SetMinimum(-1.0)
    h_ratio[hname].SetMaximum(4.0)
    
    #Switch back to first Canvas to draw the legend
    pad1[hname].cd()                       
    if signal_magnify != 1:
        h_base[myWF.sig_samplename+hname].Scale(signal_magnify)
    leg1.Draw()

    canvas[hname].SaveAs("plots/tree_" + hname + ".gif")
    #canvas[hname].SaveAs("plots/tree_" + hname + ".C")
    canvas[hname].SaveAs("plots/tree_" + hname + ".pdf")

##   Signal plots to be produced
for hname in list_histos:
    #if hname == "h_m4b" or hname == "h_m12" or  hname == "h_m34":
    if hname == "h_m4b" or hname == "h_m4b_fitted" or hname == "h_m12" or  hname == "h_m34" or hname == "h_abs_massRatio_jetpair" or hname == "h_abs_massRatio_jetpair_afterfit":
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

#canvas_mh_ma1_sig = ROOT.TCanvas("h_mh_ma1_sig","h_mh_ma1_sig")
#canvas_mh_ma1_sig.cd()
#h_mh_ma1_sig.Draw("COLZ")
#canvas_mh_ma1_sig.SaveAs("plots/tree_h_mh_ma1_signal.gif")
#canvas_mh_ma1_sig.SaveAs("plots/tree_h_mh_ma1_signal.C")
#canvas_mh_ma1_sig.SaveAs("plots/tree_h_mh_ma1_signal.pdf")

print "Number of expected events for ", luminosity_norm, " in fb-1"
print "Number of signal events passed = ", Nsig_passed
print "Number of background events passed = ", Nbkg_passed
print "Number of data events passed = ", Ndata_passed
print "Significance S/sqrt(B) = ", Nsig_passed/math.sqrt(Nbkg_passed)
print "Significance S/sqrt(B + deltaB^2) = ", Nsig_passed/(math.sqrt(Nbkg_passed) + 0.2*Nbkg_passed)
print "Significance S/sqrt(S+B) = ", Nsig_passed/math.sqrt(Nsig_passed + Nbkg_passed)
print "\nAll the intresting plots have been produced..!"
message =raw_input('\nHit the Enter Key to exit the program! ')
print(message)
print "Thank You and Good Bye.!\n\n "
exit()
