import os

###All normalizations are provided to 1fb-1 of lumi in these tables

dir_input = "crab_projects/samples/"
list_dirs = os.listdir(dir_input)

if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")

output_filename = "rootfiles/Normalizations_table.txt"

##These are in pb
def get_xsec_fromsample(samplename):
    
    if samplename == "ttbar":
        return 730.0

    if samplename == "ttbarW":
        return 0.4062
                  
    if samplename == "ttbarZ":
        return 0.5297 

    if samplename == "SingleTop_tW":
        return 35.85

    if samplename == "SingleAntiTop_tW":
        return 35.85

    if samplename == "WJetsToLNu":
        return 61526.7

    if samplename == "WJetsToQQ":
        return 95.14

    if samplename == "DY_10_50":
        return 18610.0

    if samplename == "DY_50":
        return 6104.0

    if samplename == "QCD_HT100to200":
        return 27540000.0 

    if samplename == "QCD_HT200to300_1":
        return 1717000.0

    if samplename == "QCD_HT200to300_2":
        return 1717000.0

    if samplename == "QCD_HT300to500_1":
        return 351300.0

    if samplename == "QCD_HT300to500_2":
        return 351300.0

    if samplename == "QCD_HT500to700_1":
        return 31630.0

    if samplename == "QCD_HT500to700_2":
        return 31630.0

    if samplename == "QCD_HT700to1000_1":
        return 6802.0

    if samplename == "QCD_HT700to1000_2":
        return 6802.0

    if samplename == "QCD_HT1000to1500_1":
        return 1206.0

    if samplename == "QCD_HT1000to1500_2":
        return 1206.0

    if samplename == "QCD_HT1500to2000_1":
        return 120.4 

    if samplename == "QCD_HT1500to2000_2":
        return 120.4

    if samplename == "QCD_HT2000toInf_1":
        return 25.25

    if samplename == "QCD_HT2000toInf_2":
        return 25.25

    if samplename == "QCD_MuEnriched_P20to30":
        return 2960198.40 

    if samplename == "QCD_MuEnriched_P30to50":
        return 1652471.46

    if samplename == "QCD_MuEnriched_P50to80":
        return 437504.10

    if samplename == "QCD_MuEnriched_P80to120":
        return 106033.66

    if samplename == "QCD_MuEnriched_P120to170":
        return 25190.52

    if samplename == "QCD_MuEnriched_P170to300":
        return 8654.49

    if samplename == "QCD_MuEnriched_P300to470":
        return 797.35

    if samplename == "QCD_MuEnriched_P470to600":
        return 79.03

    if samplename == "QCD_MuEnriched_P600to800":
        return 25.10

    if samplename == "QCD_MuEnriched_P800to1000":
        return 4.71

    if samplename == "QCD_MuEnriched_P1000toInf":
        return 1.62

    if samplename == "ZZ":
        return 8.16

    if samplename == "WW":
        return 51.723

    if samplename == "WZ":
        return 47.13

    if samplename == "Signal_H500_A200":
        return 1.96

    if samplename == "Signal_H700_A200":
        return 0.22

    if samplename == "Signal_H800_A300":
        return 0.10

##Now starts the program

out_file = open(output_filename,"w")

for dirname in list_dirs:
    samplename = dirname.split("crab_HAA4bAnalysis_")[1]
    print "Processing sample dir " + dirname
    crab_command = "crab report -d " + dir_input + dirname + " | grep read"
    print crab_command

    #Special cases
    #if samplename == "QCD_15_30":
    #    number_events = 38425945.*186./187.
    #    print "No. of events processed = " + str (number_events) + "\n"
    #elif samplename == "WJetsToLNu":
    #    number_events = 72207128.0
    #else :
    #    event_string = os.popen(crab_command).read()
    #    #number_events = float((event_string.split())[0])
    #    number_events = float((event_string.split())[4])
    #    print "No. of events processed = " + str (number_events) + "\n"

    event_string = os.popen(crab_command).read()
    number_events = float((event_string.split())[4])
    print "No. of events processed = " + str (number_events) + "\n"
    #xsection = float(get_xsec_fromsample(samplename))
    xsection = get_xsec_fromsample(samplename)  
    print "crsoection = ", xsection
    scale_factor = float(xsection*35.87/number_events)
    #scale_factor = float(xsection*1000./number_events)
    print "scale_factor = ", scale_factor
    write_string = samplename + " " + str(scale_factor) + "\n"
    print "Output Norm = ", write_string
    out_file.write(write_string)

print "All done!"
