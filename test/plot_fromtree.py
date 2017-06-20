import ROOT
import os
import math
import numpy as np
#import matplotlib.pyplot as plt
from ROOT import TLatex as latex

#Supress the opening of many Canvas's
ROOT.gROOT.SetBatch(True)   

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
PT1_MIN = 130.
PT2_MIN = 90.
PT3_MIN = 80.
PT4_MIN = 50.
ETA1_MAX = 5.
ETA2_MAX = 5.
ETA3_MAX = 5.
ETA4_MAX = 5.

DELTA_PHI_MIN = 0.
DELTA_ETA_MAX = 2.

PT_PAIR1_MIN = 0.
PT_PAIR2_MIN = 0.

DELTA_MASS = 0.9

JET1_BTAG = 0.8
JET2_BTAG = 0.8
JET3_BTAG = 0.8
JET4_BTAG = 0.8

#data legend name
data_legend_name = "Data"

#Normalize to this luminsity, in fb-1
#luminosity_norm = 36.46
luminosity_norm = 35.87

#Make signal histos larger
signal_magnify = 100.

output_dir = "plots"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

list_histos = ["h_jet1pt", "h_jet2pt", "h_jet3pt", "h_jet4pt", "h_delta_Phi_pair", "h_delta_Eta_pair", "h_pt_pair1", "h_pt_pair2", "h_m4b", "h_m4b_fitted", "h_jet1eta", "h_jet2eta", "h_jet3eta", "h_jet4eta", "h_jet1Btag", "h_jet2Btag", "h_jet3Btag", "h_jet4Btag", "h_m12", "h_m34","h_nPv" , "h_abs_massRatio_jetpair_beforefit", "h_abs_massRatio_jetpair_afterfit" , "h_GammaL2_X_GammaL3", "h_Gamma_triple" ] #"h_njetE" 

# Dataset names, and add new data samples as many as available
data_samples = ["BTagCSV_B", "BTagCSV_C", "BTagCSV_D", "BTagCSV_E", "BTagCSV_F", "BTagCSV_G"] 

signal_histos = ["h_m4b","h_m4b_fitted", "h_m12", "h_m34", "h_abs_massRatio_jetpair_beforefit", "h_abs_massRatio_jetpair_afterfit"] # This is just to plot the signal histograms separately.

colors_mask = [1,400,840,616,860,432,880,416,800,900,820,920]   

def select_all_but_one(cutstring):

    btag_collector = [mytree.jet1Btag, mytree.jet2Btag, mytree.jet3Btag, mytree.jet4Btag]

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

    selection_bools["h_jet1Btag"] = btag_collector[0] > JET1_BTAG
    selection_bools["h_jet2Btag"] = btag_collector[1] > JET2_BTAG
    selection_bools["h_jet3Btag"] = btag_collector[2] > JET3_BTAG
    selection_bools["h_jet4Btag"] = btag_collector[3] > JET4_BTAG

    selection_bools["h_abs_massRatio_jetpair_beforefit"] = abs_massRatio_jetpair < DELTA_MASS

    btag_select = False
    btag_counter = 0
    for btag_value in btag_collector:
        if btag_value > 0.935:
            btag_counter += 1

    result = True

    for hname in selection_bools:
        if cutstring == hname:
            continue
        else:
            result = result and selection_bools[hname]

    return result and btag_counter > 1

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

##Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##Get data file names and the corresponding rootfiles
dataName_list = myWF.get_dataSample_names()
root_file_data = myWF.get_data_root_files()

isDataAbsent = False
if not dataName_list:
    isDataAbsent = True
    print "No data available, will only use MC"

#Combining two input lists i,e sample monte carlo and data
#combined_list = dict()
combined_list = samplename_list + dataName_list

# Sorting list
combined_list.sort()  # Just to add the data entries at the top of the legend

#Store root file in a common list
#rootfiles_combined_list = dict()
rootfiles_combined_list = dict(root_file, **root_file_data)

##Get the handlers for all the histos
hs          = dict()
hs_data     = dict()
h_sum_mc    = dict()
h_sum_qcd   = dict()
h_sum_total = dict()      #total sum of mc and QCD
h_ratio     = dict()
h_base      = dict()
h_Signal    = dict()

#TLatex text
#latex.SetTextAlign(20)
#latex.SetTextSize(0.025)

for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")
    hs_data[hname] = ROOT.THStack("hs_"+ hname, "")
##Define the histos to be created

isQCDfirst = True
#for sample_name in samplename_list:   #for background only
#for sample_name in dataName_list:     #for data only
for sample_name in combined_list:      #for data + background

    if "QCD" in sample_name:
        theSampleName = "QCD_"
        if not isQCDfirst:
            continue
        isQCDfirst = False
    elif sample_name in data_samples:
        theSampleName = "Data_"
    else:
        theSampleName = sample_name

    if "Signal" in sample_name and not sample_name == myWF.sig_samplename:
        continue

    h_base[theSampleName+list_histos[0]]  = ROOT.TH1F(theSampleName+list_histos[0], "p_{T} of the 1st jet", 55, PT1_MIN, 390.)
    h_base[theSampleName+list_histos[1]]  = ROOT.TH1F(theSampleName+list_histos[1], "p_{T} of the 2nd jet", 55, PT2_MIN, 400.)
    h_base[theSampleName+list_histos[2]]  = ROOT.TH1F(theSampleName+list_histos[2], "p_{T} of the 3rd jet", 55, PT3_MIN, 400.)
    h_base[theSampleName+list_histos[3]]  = ROOT.TH1F(theSampleName+list_histos[3], "p_{T} of the 4th jet", 55, PT4_MIN, 390.)
    h_base[theSampleName+list_histos[4]]  = ROOT.TH1F(theSampleName+list_histos[4], "#Delta#phi of the two jet pairs", 40, 0., 3.14)
    h_base[theSampleName+list_histos[5]]  = ROOT.TH1F(theSampleName+list_histos[5], "#Delta#eta of the two jet pairs", 40, -5., 5.)
    h_base[theSampleName+list_histos[6]]  = ROOT.TH1F(theSampleName+list_histos[6], "p_{T} of the first jet pair", 50, 0., 500.)
    h_base[theSampleName+list_histos[7]]  = ROOT.TH1F(theSampleName+list_histos[7], "p_{T} of the second jet pair", 50, 0., 500.)
    h_base[theSampleName+list_histos[8]]  = ROOT.TH1F(theSampleName+list_histos[8], "4jets invariant mass", 55, 150., 1400.)
    h_base[theSampleName+list_histos[9]]  = ROOT.TH1F(theSampleName+list_histos[9], "Fitted 4jets invariant mass", 55, 150., 1400.)
    h_base[theSampleName+list_histos[10]] = ROOT.TH1F(theSampleName+list_histos[10], "#eta of the 1st jet", 40, -3.5, 3.5)
    h_base[theSampleName+list_histos[11]] = ROOT.TH1F(theSampleName+list_histos[11], "#eta of the 2nd jet", 40, -3.5, 3.5)
    h_base[theSampleName+list_histos[12]] = ROOT.TH1F(theSampleName+list_histos[12], "#eta of the 3rd jet", 40, -3.5, 3.5)
    h_base[theSampleName+list_histos[13]] = ROOT.TH1F(theSampleName+list_histos[13], "#eta of the 4th jet", 40, -3.5, 3.5)
    h_base[theSampleName+list_histos[14]] = ROOT.TH1F(theSampleName+list_histos[14], "B-tag of the 1st jet", 40, JET1_BTAG+0.09, 1.)
    h_base[theSampleName+list_histos[15]] = ROOT.TH1F(theSampleName+list_histos[15], "B-tag of the 2nd jet", 40, JET2_BTAG+0.04, 1.)
    h_base[theSampleName+list_histos[16]] = ROOT.TH1F(theSampleName+list_histos[16], "B-tag of the 3rd jet", 40, JET3_BTAG, 1.)
    h_base[theSampleName+list_histos[17]] = ROOT.TH1F(theSampleName+list_histos[17], "B-tag of the 4th jet", 40, JET4_BTAG, 1.)
    h_base[theSampleName+list_histos[18]] = ROOT.TH1F(theSampleName+list_histos[18], "Invariant mass m_{12} of the first jet pair", 55, 0., 650.)
    h_base[theSampleName+list_histos[19]] = ROOT.TH1F(theSampleName+list_histos[19], "Invariant mass m_{34} of the second jet pair", 55, 0., 650.)
    h_base[theSampleName+list_histos[20]] = ROOT.TH1F(theSampleName+list_histos[20], "No. of primary verticies",55 , 0., 45.)
    h_base[theSampleName+list_histos[21]] = ROOT.TH1F(theSampleName+list_histos[21], "abs([M(b1b2)-M(b3b4)]/[M(b1b2)+M(b3b4])",40 , 0., 1.1)
    h_base[theSampleName+list_histos[22]] = ROOT.TH1F(theSampleName+list_histos[22], "abs([M(b1b2)-M(b3b4)]/[M(b1b2)+M(b3b4])",40 , 0., 1.1)
    h_base[theSampleName+list_histos[23]] = ROOT.TH1F(theSampleName+list_histos[23], "#gamma_{L2} x #gamma_{L3}",50 , 0., 35.)
    h_base[theSampleName+list_histos[24]] = ROOT.TH1F(theSampleName+list_histos[24], "#gamma_{tripple}",50 , 0., 3.5)
    #h_base[theSampleName+list_histos[20]] = ROOT.TH1F(theSampleName+list_histos[20], "No. of jet entries", 10, 0., 10.)

#Defining 2D Histograms/Correlation's
h_ma1_ma2     = ROOT.TH2F("h_ma1_ma2", "m_{12} vs m_{34}", 25,   0.,  600., 25, 0., 600.)
h_mh_ma1      = ROOT.TH2F("h_mh_ma2", "m_{H} vs m_{34}"  , 25, 200., 1000., 35, 0., 600.)
h_ma1_ma2_sig = ROOT.TH2F("h_ma1_ma2_sig", "m_{12} vs m_{34}", 25,   0.,  600., 25, 0., 600.)
h_mh_ma1_sig  = ROOT.TH2F("h_mh_ma2_sig", "m_{H} vs m_{34}",   25, 200., 1000., 25, 0., 600.)

##Graphics stuff
canvas = dict()
canvas_sig = dict()
pad1 = dict()
pad2 = dict()

for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname,200,106,600,600)
    canvas_sig[hname] = ROOT.TCanvas(hname+"_sig",hname+"_sig",200,106,600,600)

    pad1[hname] = ROOT.TPad(hname, hname,0,0.2097902,0.9966443,0.9667832)   #changes from here
    pad2[hname] = ROOT.TPad(hname, hname, 0,0.005244755,0.9966443,0.2465035)    
    
leg1 = ROOT.TLegend(0.6868687,0.6120093,0.9511784,0.9491917)
leg1.SetHeader(" ")
leg1.SetFillColor(0)
leg1.SetBorderSize(0)
leg1.SetLineColor(1)
leg1.SetLineStyle(1)
leg1.SetLineWidth(1)
leg1.SetFillStyle(1001)

Nsig_passed = 0.
Nbkg_passed = 0.
Ndata_passed = 0.

##Loop on samples, and then on events, and merge QCD stuff
idx_sample = 0
isFirstQCDlegend = True

#for name_sample in samplename_list: # for only background
#for  name_sample in dataName_list:  # for only data
for  name_sample in combined_list:    # for data + background

    if "Signal" in name_sample and not name_sample == myWF.sig_samplename:
        continue

    theSampleName = name_sample

    if "QCD" in name_sample:
        QCDflag = True
        theSampleName = "QCD_"
    else:
        QCDflag = False

    if name_sample in data_samples:
        Dataflag = True
        theSampleName = "Data_"
    else:
        Dataflag = False
     
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

        if Dataflag:
            Event_Weight = 1.0                     #to be used in the future for corrections
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
        if combination_flag_fit==1:
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

        pt_pair1 = p_pair1.Pt()
        pt_pair2 = p_pair2.Pt()

        totaljets_4mom = p_pair1 + p_pair2
        m4b = totaljets_4mom.M()

        pt_pair1_fit = p_pair1_fit.Pt()
        pt_pair2_fit = p_pair2_fit.Pt()

        totaljets_4mom_fitted = p_pair1_fit + p_pair2_fit
        m4b_fitted = totaljets_4mom_fitted.M()

        #abs(mass_diff_pair1/mass_add_pair2) variable before fit
        diff_p_pair = p_pair1 - p_pair2
        add_p_pair  = p_pair1 + p_pair2
        mass_diff_pair = diff_p_pair.M()
        mass_add_pair  = add_p_pair.M()
        if not mass_add_pair==0:
            abs_massRatio_jetpair = abs(mass_diff_pair/mass_add_pair)
    
        #abs(mass_diff_pair1/mass_add_pair2) variable after fit
        diff_p_pair_fit = p_pair1_fit - p_pair2_fit
        add_p_pair_fit = p_pair1_fit + p_pair2_fit
        mass_diff_pair_fit = diff_p_pair_fit.M()
        mass_add_pair_fit = add_p_pair_fit.M()

        #Define Gamma variables using energy and mass
        #double jet boast
        if not ((jet1_4mom + jet2_4mom).M()) == 0:
            Gamma12 = (jet1_4mom.Energy() + jet2_4mom.Energy())/((jet1_4mom + jet2_4mom).M())
        if not ((jet1_4mom + jet3_4mom).M()) == 0:
            Gamma13 = (jet1_4mom.Energy() + jet3_4mom.Energy())/((jet1_4mom + jet3_4mom).M())
        if not ((jet1_4mom + jet4_4mom).M()) == 0:
            Gamma14 = (jet1_4mom.Energy() + jet4_4mom.Energy())/((jet1_4mom + jet4_4mom).M())
        if not ((jet2_4mom + jet3_4mom).M()) == 0:
            Gamma23 = (jet2_4mom.Energy() + jet3_4mom.Energy())/((jet2_4mom + jet3_4mom).M())
        if not ((jet2_4mom + jet4_4mom).M()) == 0:
            Gamma24 = (jet2_4mom.Energy() + jet4_4mom.Energy())/((jet2_4mom + jet4_4mom).M())
        if not ((jet3_4mom + jet4_4mom).M()) == 0:
            Gamma34 = (jet3_4mom.Energy() + jet4_4mom.Energy())/((jet3_4mom + jet4_4mom).M())

        diff_Gamma1 = abs(Gamma12 - Gamma34)
        diff_Gamma2 = abs(Gamma13 - Gamma24)
        diff_Gamma3 = abs(Gamma14 - Gamma23)

        #Ordering as per minimum difference (Increasing order)
        # diff_Gamma1 < (diff_Gamma2 < diff_Gamma3 or diff_Gamma2 > diff_Gamma3)
        if diff_Gamma1 < diff_Gamma2 and diff_Gamma1 < diff_Gamma3:
            if diff_Gamma2 < diff_Gamma3:
                if Gamma12 > Gamma34:
                     Gamma_L1 = Gamma12
                     Gamma_S1 = Gamma34
                else:
                     Gamma_L1 = Gamma34
                     Gamma_S1 = Gamma12

                if Gamma13 > Gamma24:
                     Gamma_L2 = Gamma13
                     Gamma_S2 = Gamma24
                else:
                     Gamma_L2 = Gamma24
                     Gamma_S2 = Gamma13

                if Gamma14 > Gamma23:
                    Gamma_L3 = Gamma14
                    Gamma_S3 = Gamma23
                else:
                    Gamma_L3 = Gamma23
                    Gamma_S3 = Gamma14

            elif diff_Gamma3 < diff_Gamma2:
                if Gamma12 > Gamma34:
                     Gamma_L1 = Gamma12
                     Gamma_S1 = Gamma34
                else:
                     Gamma_L1 = Gamma34
                     Gamma_S1 = Gamma12

                if Gamma14 > Gamma23:
                    Gamma_L2 = Gamma14
                    Gamma_S2 = Gamma23
                else:
                    Gamma_L2 = Gamma23
                    Gamma_S2 = Gamma14

                if Gamma13 > Gamma24:
                     Gamma_L3 = Gamma13
                     Gamma_S3 = Gamma24
                else:
                     Gamma_L3 = Gamma24
                     Gamma_S3 = Gamma13

        # diff_Gamma2 < (diff_Gamma1 < diff_Gamma3 or diff_Gamma1> diff_Gamma3)
        if diff_Gamma2 < diff_Gamma1 and diff_Gamma2 < diff_Gamma3:
            if diff_Gamma1 < diff_Gamma3:

                if Gamma13 > Gamma24:
                     Gamma_L1 = Gamma13
                     Gamma_S1 = Gamma24
                else:
                     Gamma_L1 = Gamma24
                     Gamma_S1 = Gamma13

                if Gamma12 > Gamma34:
                     Gamma_L2 = Gamma12
                     Gamma_S2 = Gamma34
                else:
                     Gamma_L2 = Gamma34
                     Gamma_S2 = Gamma12

                if Gamma14 > Gamma23:
                     Gamma_L3 = Gamma14
                     Gamma_S3 = Gamma23
                else:
                     Gamma_L3 = Gamma23
                     Gamma_S3 = Gamma14

            elif diff_Gamma3 < diff_Gamma1:

                if Gamma13 > Gamma24:
                     Gamma_L1 = Gamma13
                     Gamma_S1 = Gamma24
                else:
                     Gamma_L1 = Gamma24
                     Gamma_S1 = Gamma13

                if Gamma14 > Gamma23:
                     Gamma_L2 = Gamma14
                     Gamma_S2 = Gamma23
                else:
                     Gamma_L2 = Gamma23
                     Gamma_S2 = Gamma14

                if Gamma12 > Gamma34:
                     Gamma_L3 = Gamma12
                     Gamma_S3 = Gamma34
                else:
                     Gamma_L3 = Gamma34
                     Gamma_S3 = Gamma12

        # diff_Gamma3 < (diff_Gamma1 < diff_Gamma2 or diff_Gamma1 > diff_Gamma2)
        if diff_Gamma3 < diff_Gamma1 and diff_Gamma3 < diff_Gamma2:
            if diff_Gamma1 < diff_Gamma2:
                if Gamma14 > Gamma23:
                     Gamma_L1 = Gamma14
                     Gamma_S1 = Gamma23
                else:
                     Gamma_L1 = Gamma23
                     Gamma_S1 = Gamma14

                if Gamma12 > Gamma34:
                     Gamma_L2 = Gamma12
                     Gamma_S2 = Gamma34
                else:
                     Gamma_L2 = Gamma34
                     Gamma_S2 = Gamma12

                if Gamma13 > Gamma24:
                     Gamma_L3 = Gamma13
                     Gamma_S3 = Gamma24
                else:
                     Gamma_L3 = Gamma24
                     Gamma_S3 = Gamma13

            elif diff_Gamma2 < diff_Gamma1:
                if Gamma14 > Gamma23:
                     Gamma_L1 = Gamma14
                     Gamma_S1 = Gamma23
                else:
                     Gamma_L1 = Gamma23
                     Gamma_S1 = Gamma14

                if Gamma13 > Gamma24:
                     Gamma_L2 = Gamma13
                     Gamma_S2 = Gamma24
                else:
                     Gamma_L2 = Gamma24
                     Gamma_S2 = Gamma13

                if Gamma12 > Gamma34:
                     Gamma_L3 = Gamma12
                     Gamma_S3 = Gamma34
                else:
                     Gamma_L3 = Gamma34
                     Gamma_S3 = Gamma12
        
        #Parameter to controll QCD multijet
        Gamma_L2_X_GammaL3 = Gamma_L2 * Gamma_L3

        #Triple Jet boast
        if not ((jet1_4mom + jet2_4mom + jet3_4mom).M()) == 0:
            Gamma123 = (jet1_4mom.Energy() + jet2_4mom.Energy() + jet3_4mom.Energy())/((jet1_4mom + jet2_4mom + jet3_4mom).M())
        if not ((jet1_4mom + jet2_4mom + jet4_4mom).M()) == 0:
            Gamma124 = (jet1_4mom.Energy() + jet2_4mom.Energy() + jet4_4mom.Energy())/((jet1_4mom + jet2_4mom + jet4_4mom).M())
        if not ((jet1_4mom + jet3_4mom + jet4_4mom).M()) == 0:
            Gamma134 = (jet1_4mom.Energy() + jet3_4mom.Energy() + jet4_4mom.Energy())/((jet1_4mom + jet3_4mom + jet4_4mom).M())
        if not ((jet2_4mom + jet3_4mom + jet4_4mom).M()) == 0:
            Gamma234 = (jet2_4mom.Energy() + jet3_4mom.Energy() + jet4_4mom.Energy())/(((jet2_4mom + jet3_4mom + jet4_4mom).M()))

        #find the maximal of the above boast
        if Gamma123 > Gamma124 and Gamma123 > Gamma134 and Gamma123 > Gamma234:
            Gamma_tripple = Gamma123
        if Gamma124 > Gamma123 and Gamma124 > Gamma134 and Gamma124 > Gamma234:
            Gamma_tripple = Gamma124
        if Gamma134 > Gamma123 and Gamma134 > Gamma124 and Gamma134 > Gamma234:
            Gamma_tripple = Gamma134
        if Gamma234 > Gamma123 and Gamma234 > Gamma124 and Gamma234 > Gamma134:
            Gamma_tripple = Gamma234


        if not mass_add_pair_fit == 0:
            abs_massRatio_jetpair_afterfit = abs(mass_diff_pair_fit/mass_add_pair_fit)

        if select_all_but_one("h_jet1pt"):
            h_base[theSampleName+"h_jet1pt"].Fill(jet1pt,Event_Weight)
        if select_all_but_one("h_jet2pt"):
            h_base[theSampleName+"h_jet2pt"].Fill(jet2pt,Event_Weight)
        if select_all_but_one("h_jet3pt"):
            h_base[theSampleName+"h_jet3pt"].Fill(jet3pt,Event_Weight)
        if select_all_but_one("h_jet4pt"):
            h_base[theSampleName+"h_jet4pt"].Fill(jet4pt,Event_Weight)
        if select_all_but_one("h_delta_Phi_pair"):
            h_base[theSampleName+"h_delta_Phi_pair"].Fill(delta_Phi,Event_Weight)
        if select_all_but_one("h_delta_Eta_pair"):
            h_base[theSampleName+"h_delta_Eta_pair"].Fill(delta_Eta,Event_Weight)
        if select_all_but_one("h_pt_pair1"):
            h_base[theSampleName+"h_pt_pair1"].Fill(pt_pair1,Event_Weight)
        if select_all_but_one("h_pt_pair2"):
            h_base[theSampleName+"h_pt_pair2"].Fill(pt_pair2,Event_Weight)
        if select_all_but_one("h_m4b"):
            h_base[theSampleName+"h_m4b"].Fill(m4b,Event_Weight)
        if select_all_but_one("h_m4b_fitted"):
            h_base[theSampleName+"h_m4b_fitted"].Fill(m4b_fitted,Event_Weight)
        if select_all_but_one("h_jet1eta"):
            h_base[theSampleName+"h_jet1eta"].Fill(jet1eta,Event_Weight)
        if select_all_but_one("h_jet2eta"):
            h_base[theSampleName+"h_jet2eta"].Fill(jet2eta,Event_Weight)
        if select_all_but_one("h_jet3eta"):
            h_base[theSampleName+"h_jet3eta"].Fill(jet3eta,Event_Weight)
        if select_all_but_one("h_jet4eta"):
            h_base[theSampleName+"h_jet4eta"].Fill(jet4eta,Event_Weight)
        if select_all_but_one("h_jet1Btag"):
            h_base[theSampleName+"h_jet1Btag"].Fill(mytree.jet1Btag,Event_Weight)
        if select_all_but_one("h_jet2Btag"):
            h_base[theSampleName+"h_jet2Btag"].Fill(mytree.jet2Btag,Event_Weight)
        if select_all_but_one("h_jet3Btag"):
            h_base[theSampleName+"h_jet3Btag"].Fill(mytree.jet3Btag,Event_Weight)
        if select_all_but_one("h_jet4Btag"):
            h_base[theSampleName+"h_jet4Btag"].Fill(mytree.jet4Btag,Event_Weight)
        if select_all_but_one("h_m12"): 
            h_base[theSampleName+"h_m12"].Fill(p_pair1.M(),Event_Weight)
        if select_all_but_one("h_m34"):
            h_base[theSampleName+"h_m34"].Fill(p_pair2.M(),Event_Weight)
        if select_all_but_one("h_nPv"):
            h_base[theSampleName+"h_nPv"].Fill(mytree.N_nPv,Event_Weight)
        if select_all_but_one("h_abs_massRatio_jetpair_beforefit"):
            h_base[theSampleName+"h_abs_massRatio_jetpair_beforefit"].Fill(abs_massRatio_jetpair,Event_Weight)
        if select_all_but_one("h_abs_massRatio_jetpair_afterfit"):
            h_base[theSampleName+"h_abs_massRatio_jetpair_afterfit"].Fill(abs_massRatio_jetpair_afterfit, Event_Weight)
        if select_all_but_one("h_Gamma_L2_X_GammaL3"):
            h_base[theSampleName+"h_GammaL2_X_GammaL3"].Fill(Gamma_L2_X_GammaL3,Event_Weight)
        if select_all_but_one("h_Gamma_triple"):
            h_base[theSampleName+"h_Gamma_triple"].Fill(Gamma_tripple,Event_Weight)

        ##Now the 2D plots
        if name_sample == myWF.sig_samplename:
            h_ma1_ma2_sig.Fill(p_pair1.M(),p_pair2.M(),Event_Weight)
        elif not Dataflag:
            h_ma1_ma2.Fill(p_pair1.M(),p_pair2.M(),Event_Weight)

        if name_sample == myWF.sig_samplename:
            h_mh_ma1_sig.Fill(m4b,p_pair1.M(),Event_Weight)   
        elif not Dataflag:
            h_mh_ma1.Fill(m4b,p_pair1.M(),Event_Weight)   

        #Count the events
        if select_all_but_one("actually all"):
            if name_sample == myWF.sig_samplename:
                Nsig_passed += norm_factor
            elif Dataflag:
                Ndata_passed += norm_factor
            else:
                Nbkg_passed += norm_factor

    for idx_histo,hname in enumerate(list_histos):

        if QCDflag:
            h_base[theSampleName+hname].SetFillColor(2)
        elif Dataflag:
            h_base[theSampleName+hname].SetMarkerStyle(20)
            h_base[theSampleName+hname].SetMarkerColor(1)
            h_base[theSampleName+hname].SetMarkerSize(1)
        elif name_sample == myWF.sig_samplename:
            h_base[theSampleName+hname].SetLineStyle(2)   #dashed
            h_base[theSampleName+hname].SetLineColor(4)   #blue
            h_base[theSampleName+hname].SetLineWidth(4)   #kind of thick
        else:
            h_base[theSampleName+hname].SetFillColor(colors_mask[idx_sample])

        if idx_histo == 0:
            if Dataflag:
                leg1.AddEntry(h_base[theSampleName+hname],data_legend_name,"elp")               
            elif QCDflag and isFirstQCDlegend:
                 leg1.AddEntry(h_base[theSampleName+hname],"QCD","f")
                 isFirstQCDlegend = False
            elif name_sample == myWF.sig_samplename:
                 sample_legend_name = "100 x " + name_sample
                 leg1.AddEntry(h_base[name_sample+hname], sample_legend_name,"f")  #To comment when signal is has to be excluded.
            elif not QCDflag:
                 leg1.AddEntry(h_base[theSampleName+hname],theSampleName,"f")

        #elif name_sample == myWF.sig_samplename:    #To uncomment it when the signal needs to be excluded
        #     continue
        if not QCDflag and not Dataflag:
            hs[hname].Add(h_base[theSampleName+hname])
            h_sum_mc[hname] = ROOT.TH1F(h_base[theSampleName+hname])
            h_sum_mc[hname].Add(h_base[theSampleName+hname])        

    if not QCDflag and not Dataflag and not name_sample == myWF.sig_samplename:
        idx_sample += 1

print "Finished runnning over samples!"

for idx_histo,hname in enumerate(list_histos):
    hs[hname].Add(h_base["QCD_"+hname])
    h_sum_qcd[hname]=ROOT.TH1F(h_base["QCD_"+hname]) 
    h_sum_qcd[hname].Add(h_base["QCD_"+hname])

for hname in list_histos: 
    canvas[hname].cd()
    #canvas[hname].SetLogy(1)

    #Draw two pads pad1 and pad2
    pad1[hname].Draw()                      
    pad1[hname].SetFillColor(0)
    pad1[hname].SetFrameFillColor(0)
    pad1[hname].SetBorderSize(2)
    pad1[hname].SetRightMargin(0.04545)
    pad1[hname].SetTopMargin(0.0369)
    pad1[hname].SetFrameBorderMode(0)
   
    pad2[hname].Draw()
    pad2[hname].SetBorderMode(0)
    pad2[hname].SetBorderSize(2)
    pad2[hname].SetGridy()
    pad2[hname].SetRightMargin(0.04545)
    pad2[hname].SetTopMargin(0.01493)
    pad2[hname].SetBottomMargin(0.269)
    pad2[hname].SetFrameBorderMode(0)
    
    pad2[hname].SetFrameFillColor(0)
    pad2[hname].SetFillColor(0)
    
    #Switching to first pad in the same canvas
    pad1[hname].cd()                          
    hs[hname].Draw("histo")

    #Draw the data separately
    if not isDataAbsent:
        h_base["Data_"+hname].Draw("same e1p") 
        hmax_data = h_base["Data_"+hname].GetMaximum()
        hmax = hs[hname].GetMaximum()

        if hmax < hmax_data:
            hmax = hmax_data

        hs[hname].SetMaximum(200 + hmax)
        h_base["Data_"+hname].SetMaximum(200 + hmax)
      
    #Grapic names
    hs[hname].SetTitle(" ")
    hs[hname].GetXaxis().SetTickLength(0.03) 
    hs[hname].GetXaxis().SetLabelOffset(0.006)
    hs[hname].GetYaxis().SetLabelOffset(0.007)
    hs[hname].GetYaxis().SetTitleOffset(1.077) #####
    hs[hname].GetYaxis().SetLabelSize(0.048)
    hs[hname].GetYaxis().SetNdivisions(506)
    hs[hname].GetYaxis().SetTitle("Events")
    hs[hname].GetYaxis().SetTitleSize(0.045)
    hs[hname].GetXaxis().SetLabelFont(0)
    hs[hname].GetYaxis().CenterTitle()

    if not isDataAbsent:
        h_base["Data_"+hname].GetYaxis().SetLabelSize(0.048)
        h_base["Data_"+hname].GetYaxis().SetNdivisions(506)

    #Switching to second pad on the same canvas
    pad2[hname].cd() 
    ROOT.gStyle.SetOptStat(0)

    if not isDataAbsent:
        h_sum_total[hname]=ROOT.TH1F(h_sum_mc[hname])
        h_sum_total[hname].Add(h_sum_qcd[hname])

        #Take the ratio between data and MC and draw it
        h_ratio[hname]= ROOT.TH1F(h_base["Data_"+hname])
        h_ratio[hname].Divide(h_sum_total[hname])                     
        h_ratio[hname].Draw("elp")

        #Graphic names 
        h_ratio[hname].GetXaxis().SetTitle(h_base[name_sample+hname].GetTitle())
        h_ratio[hname].SetTitle(" ")
        h_ratio[hname].GetXaxis().SetTickLength(0.04)    
        h_ratio[hname].GetXaxis().SetLabelFont(0)     
        h_ratio[hname].GetYaxis().SetNdivisions(506)     
        h_ratio[hname].GetYaxis().SetLabelSize(0.12)
        #h_ratio[hname].GetXaxis().SetLabelSize(0.50)   
        h_ratio[hname].GetXaxis().SetLabelSize(0.125)     
        h_ratio[hname].GetXaxis().SetTickLength(0.11)
        h_ratio[hname].GetYaxis().SetTitle("Data/MC")
        h_ratio[hname].GetYaxis().SetTitleSize(0.125)
        h_ratio[hname].GetYaxis().SetTitleOffset(0.30)
        h_ratio[hname].GetXaxis().SetTitleOffset(0.80)
        h_ratio[hname].GetXaxis().SetTitleSize(0.15)
        ROOT.gPad.SetGridy(1)
        h_ratio[hname].GetXaxis().CenterTitle()
        h_ratio[hname].GetYaxis().CenterTitle()
        h_ratio[hname].SetMinimum(0.0)
        h_ratio[hname].SetMaximum(2.0)
    
    #Switch back to first Canvas to draw the legend
    pad1[hname].cd()                       
    if signal_magnify != 1:
        h_base[myWF.sig_samplename+hname].Scale(signal_magnify)
    leg1.Draw()

    canvas[hname].SaveAs("plots/tree_" + hname + ".gif")
    #canvas[hname].SaveAs("plots/tree_" + hname + ".C")
    canvas[hname].SaveAs("plots/tree_" + hname + ".pdf")

#Signal only plots to be produced
for hname in list_histos:
    if hname in signal_histos: 
        canvas_sig[hname].cd()

        ##Remove some style fancyness
        h_base[myWF.sig_samplename+hname].SetFillColor(0)
        h_base[myWF.sig_samplename+hname].SetLineStyle(0)

        ROOT.gStyle.SetOptStat(1)
        h_base[myWF.sig_samplename+hname].Draw("histo")
        ##################Naming
        ##h_base[myWF.sig_samplename+hname].SetTitle("Signal Plots")
        h_base[myWF.sig_samplename+hname].GetYaxis().SetTitle("Events")
        h_base[myWF.sig_samplename+hname].GetXaxis().SetTitle(h_base[myWF.sig_samplename+hname].GetTitle())
        h_base[myWF.sig_samplename+hname].SetTitle("Signal")
        h_base[myWF.sig_samplename+hname].GetXaxis().CenterTitle()
        h_base[myWF.sig_samplename+hname].GetYaxis().CenterTitle()
        h_base[myWF.sig_samplename+hname].GetXaxis().SetTickLength(0.03) 
        h_base[myWF.sig_samplename+hname].GetXaxis().SetLabelOffset(0.006)
        h_base[myWF.sig_samplename+hname].GetYaxis().SetLabelOffset(0.007)
        h_base[myWF.sig_samplename+hname].GetYaxis().SetLabelSize(0.048)
        h_base[myWF.sig_samplename+hname].GetYaxis().SetNdivisions(506)
        h_base[myWF.sig_samplename+hname].GetYaxis().SetTitleSize(0.045)
      
        #####################
        canvas_sig[hname].SaveAs("plots/tree_" + hname + "_signal.gif")
        #canvas_sig[hname].SaveAs("plots/tree_" + hname + "_signal.C")
        canvas_sig[hname].SaveAs("plots/tree_" + hname + "_signal.pdf")

##Now the 2D plots
ROOT.gStyle.SetOptStat(0)

canvas_ma1_ma2 = ROOT.TCanvas("h_ma1_ma2","h_ma1_ma2")
canvas_ma1_ma2.cd()
h_ma1_ma2.Draw("COLZ")
canvas_ma1_ma2.SaveAs("plots/tree_h_ma1_ma2.gif")
canvas_ma1_ma2.SaveAs("plots/tree_h_ma1_ma2.pdf")

canvas_ma1_ma2_sig = ROOT.TCanvas("h_ma1_ma2_sig","h_ma1_ma2_sig")
canvas_ma1_ma2_sig.cd()
h_ma1_ma2_sig.Draw("COLZ")
canvas_ma1_ma2_sig.SaveAs("plots/tree_h_ma1_ma2_signal.gif")
canvas_ma1_ma2_sig.SaveAs("plots/tree_h_ma1_ma2_signal.pfd")

canvas_mh_ma1 = ROOT.TCanvas("h_mh_ma1","h_mh_ma1")
canvas_mh_ma1.cd()
h_mh_ma1.Draw("COLZ")
canvas_mh_ma1.SaveAs("plots/tree_h_mh_ma1.gif")
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
