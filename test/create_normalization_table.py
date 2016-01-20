import ROOT
import os

###All normalizations are provided to 1fb-1 of lumi in these tables

dir_input = "crab_projects/backgrounds/"
list_dirs = os.listdir(dir_input)

output_filename = "Background_normalizations.txt"

##These are in pb
def get_xsec_fromsample(samplename):
    
    if samplename == "ttbar":
        return 670.3

    if samplename == "DY_5_50":
        return 71600.0

    if samplename == "QCD_120_170":
        return 471100.0

    if samplename == "QCD_170_300":
        return 117276.0

    if samplename == "QCD_300_470":
        return 7823.0

    if samplename == "QCD_470_600":
        return 648.2

    if samplename == "QCD_600_800":
        return 186.9

    if samplename == "QCD_800_1000":
        return 32.293

    if samplename == "QCD_1000_1400":
        return 9.4183

    if samplename == "QCD_1400_1800":
        return 0.84265

    if samplename == "QCD_1800_2400":
        return 0.114943

    if samplename == "SingleTop_tW":
        return 38.09

    if samplename == "SingleAntiTop_tW":
        return 38.09

    if samplename == "ZZ":
        return 10.32

    if samplename == "WW":
        return 63.21

##Now starts the program

out_file = open(output_filename,"w")


for dirname in list_dirs:

    samplename = dirname.split("crab_HAA4bAnalysis_")[1]
    print "Processing sample dir " + dirname
    crab_command = "crab report -d " + dir_input + dirname + "| grep read"

    number_events = -1.
    if(samplename == "ttbar"):
        number_events = 42784971.0
    else :
        event_string = os.popen(crab_command).read()
        number_events = float((event_string.split())[0])

    xsection = float(get_xsec_fromsample(samplename))

    scale_factor = xsection*1000./number_events

    write_string = samplename + " " + str(scale_factor) + "\n"
    out_file.write(write_string)

out_file.close()
