import ROOT
import os

dir_input = "crab_projects/samples/"
dir_output = "rootfiles/"

list_dirs = os.listdir(dir_input)

for dirname in list_dirs:

    print "Processing sample dir " + dirname
    crab_command = "crab getoutput -d " + dir_input + dirname
    #os.system(crab_command)

    samplename = dirname.split("crab_HAA4bAnalysis_")

    hadd_command = "hadd " + dir_output + "/HAA4bAnalysis_" + samplename[1] + ".root " + dir_input + dirname + "/results/*.root"
    os.system(hadd_command)

print "All done!"
