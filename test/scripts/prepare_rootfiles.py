import ROOT
import os

dir_input = "../crab_config/crab_projects/samples/"
dir_output_bkg = "/rootfiles/backgrounds/"
dir_output_sig = "/rootfiles/signals/"

list_dirs = os.listdir(dir_input)

if not os.path.exists(dir_output_bkg):
    os.makedirs(dir_output_bkg)

if not os.path.exists(dir_output_sig):
    os.makedirs(dir_output_sig)

for dirname in list_dirs:

    print "Processing sample dir " + dirname
    crab_command = "crab getoutput -d " + dir_input + dirname
    os.system(crab_command)

    samplename = dirname.split("crab_HAA4bAnalysis_")

    if "Signal" in dirname:
        hadd_command = "hadd " + dir_output_sig + "/HAA4bAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    else:
        hadd_command = "hadd " + dir_output_bkg + "/HAA4bAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"

    os.system(hadd_command)

print "All done!"
