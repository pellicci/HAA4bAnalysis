import ROOT
import math
import numpy as np

from Workflow_Handler import Workflow_Handler
myWF = Workflow_Handler("Signal_H500_A200")

##Global constants
PT1_MIN = 50.
DELTA_PHI_MIN = 0.

list_histos = ["h_jet1pt"]

def select_all_but_one(cutstring):

    selection_bools = dict()
    selection_bools["h_jet1pt"] = mytree.jet1pt > PT1_MIN
    selection_bools["h_delta_Phi_pair"] = mytree.delta_Phi_pair > DELTA_PHI_MIN

    result = True

    for hname in selection_bools:
        if cutstring == hname:
            continue
        else result = result or selection_bools[hname]

    return result

##Here starts the program
Norm_Map = myWF.get_normalizations_map()

##Get the handlers for all the histos
hs = dict()
for hname in list_histos:
    hs[hname] = ROOT.THStack("hs_" + hname,"")

##Define the histos to be created
h_base = dict()
h_base[list_histos[0]]  = ROOT.TH1F(list_histos[0], "p_{t} of 1st jet", 200, 0., 500.)

h_QCD = dict()
for hname in list_histos:
    h_QCD[hname]  = h_base[hname].Clone()

##Get the files and the names of the samples
samplename_list = myWF.get_samples_names()
root_file = myWF.get_root_files()

##Graphics stuff
canvas = dict()
for hname in list_histos:
    canvas[hname] = ROOT.TCanvas(hname,hname)
leg1 = ROOT.TLegend(0.6,0.6,0.9,0.9)
leg1.SetHeader("Samples considered")

##First, loop on the QCD samples and merge them
for name_sample in samplename_list:
    if not "QCD" in name_sample:
        continue

    norm_factor = Norm_Map[name_sample]
    mytree = root_file[name_sample].Get("HZZ4bAnalysis/mytree")

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break

        nb = mytree.GetEntry( jentry )
        if nb <= 0:
            continue

        if select_all_but_one("h_jet1pt"):
            h_QCD["h_jet1pt"].Fill(mytree.jet1pt,norm_factor)

for idx_histo,hname in enumerate(list_histos):
    hQCD[hname].SetFillColor(2) #Red
    if idx_hname == 0:
        leg1.AddEntry(h_QCD[hname],"QCD","f")
    hs[hname].Add(h_QCD[hname])

##Loop on all the other samples
for idx_sample,name_sample in enumerate(samplename_list):
    if "QCD" in name_sample:
        continue

    norm_factor = Norm_Map[name_sample]
    mytree = root_file[name_sample].Get("HZZ4bAnalysis/mytree")

    h_tmp = dict()
    for hname in list_histos:
        h_tmp[hname]  = h_base[hname].Clone()

    for jentry in xrange(mytree.GetEntriesFast()):
        ientry = mytree.LoadTree( jentry )
        if ientry < 0:
            break

        nb = mytree.GetEntry( jentry )
        if nb <= 0:
            continue

        if select_all_but_one("h_jet1pt"):
            h_tmp["h_jet1pt"].Fill(mytree.jet1pt,norm_factor)

    for hname in list_histos:
        h_tmp[hname].SetFillColor(idx_sample+3)
        if name_sample == myWF.sig_samplename:
                h_tmp[hname].SetLineStyle(2)   #dashed
                h_tmp[hname].SetLineColor(4)   #blue
                h_tmp[hname].SetLineWidth(4)   #kind of thik
        leg1.AddEntry(h_tmp[hname],name_sample,"f")
        hs[hname].Add(h_tmp[hname])

for hname in list_histos:
    canvas[hname].cd()
    hs[hname].Draw()
    leg1.Draw()
    canvas[hname].SaveAs("plots/tree_" + hname + ".gif")
