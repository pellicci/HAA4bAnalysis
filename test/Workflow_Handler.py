import ROOT
import os

class Workflow_Handler:

    def __init__(self,signalname, subprocess="//"):

        ##Where the root files are
        self.subprocess = subprocess

        self.norm_filename = "rootfiles/" + self.subprocess + "Normalizations_table.txt"
        self.dir_back_input = "rootfiles/" + self.subprocess + "backgrounds/"

        self.sig_samplename = signalname
        self.sig_filename = "rootfiles/" + self.subprocess + "signals/" + "HAA4bAnalysis_" + self.sig_samplename + ".root"

    def get_normalizations_map(self):
        in_file = open(self.norm_filename,"r")
        norm_map = dict()

        for line in in_file:
            data_norm = line.split()
            norm_map[data_norm[0]] = float(data_norm[1])

        return norm_map

    def get_samples_names(self, Add_Signal=True):
        list_dirs_bkg = os.listdir(self.dir_back_input)
        samplename_list = []

        for dirname in list_dirs_bkg:
            tmp_samplename = dirname.split("HAA4bAnalysis_")[1]
            tmp_samplename2 = tmp_samplename.replace(".root","")
            samplename_list.append(tmp_samplename2)

        if Add_Signal:
            samplename_list.append(self.sig_samplename)
        return samplename_list

    def get_root_files(self,Add_Signal=True):
        list_dirs_bkg = os.listdir(self.dir_back_input)
        root_file = dict()

        for dirname in list_dirs_bkg:
            tmp_samplename = dirname.split("HAA4bAnalysis_")[1]
            tmp_samplename2 = tmp_samplename.replace(".root","")
            root_file[tmp_samplename2] = ROOT.TFile(self.dir_back_input + dirname)

        if Add_Signal:
            root_file[self.sig_samplename] = ROOT.TFile(self.sig_filename)

        return root_file

    def get_best_combination(self, m1, m2, m3, m4):

        diff_m_12_34 = abs( (m1+m2).M() - (m3+m4).M() );
        diff_m_13_24 = abs( (m1+m3).M() - (m2+m4).M() );
        diff_m_14_23 = abs( (m1+m4).M() - (m2+m3).M() );

        diff_vector = [diff_m_12_34, diff_m_13_24, diff_m_14_23]
        diff_vector.sort()

        if diff_vector[0] == diff_m_12_34:
            return 1
        elif diff_vector[0] == diff_m_13_24:
            return 2
        elif diff_vector[0] == diff_m_14_23:
            return 3



