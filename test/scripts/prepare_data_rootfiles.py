import ROOT
import os

dir_input = "../crab_config/crab_projects/data/"
dir_output = "../rootfiles/data/"

list_dirs = os.listdir(dir_input)

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

for dirname in list_dirs:

    print "Processing sample dir " + dirname
    crab_command = "crab getoutput -d " + dir_input + dirname
    os.system(crab_command)

    samplename = dirname.split("crab_HAA4bAnalysis_")

    hadd_command = "hadd " + dir_output + "/HAA4bAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"

    os.system(hadd_command)

print "All done!"

