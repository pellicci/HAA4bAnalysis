import os

###All normalizations are provided to 1fb-1 of lumi in these tables

dir_input = "crab_projects/samples/"
list_dirs = os.listdir(dir_input)

# Just to write 1.0 for data to the output file for normalization
dir_input_data = "crab_projects/data/"
list_dirs_data = os.listdir(dir_input_data)

if not os.path.exists("rootfiles"):
    os.makedirs("rootfiles")

output_filename = "rootfiles/Normalizations_table.txt"

##These are in pb
def get_xsec_fromsample(samplename):
    
    if samplename == "DY_5_50":
        return 7160.0

    if samplename == "DY_50":
        return 4895.0

    if samplename == "DY_100_200":
        return 226.0

    if samplename == "DY_200_400":
        return 7.67

    if samplename == "ttbar":
        return 831.76

    if samplename == "SingleTop_tW":
        return 35.6

    if samplename == "SingleAntiTop_tW":
        return 35.6


    if samplename == "WJetsToLNu":
        return 61526.7

# Second Addition
    if samplename == "ttbarW":
        return 0.27

    if samplename == "ttbarZ":
        return 0.37

    if samplename == "QCD_HT100to200":
        return 27990000.0

    if samplename == "QCD_HT200to300":
        return 1712000.0

    if samplename == "QCD_HT300to500":
        return 347700.0

    if samplename == "QCD_HT500to700":
        return 32100.0

    if samplename == "QCD_HT700to1000":
        return 6831.0

    if samplename == "QCD_HT1000to1500":
        return 1207.0

    if samplename == "QCD_HT1500to2000":
        return 119.9

    if samplename == "QCD_HT2000toInf":
        return 25.24

# Third Addition

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
        return 23.50

    if samplename == "WZ":
        return 2.48 

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
    scale_factor = float(xsection*1000./number_events)
    write_string = samplename + " " + str(scale_factor) + "\n"
    out_file.write(write_string)

#====================================================================================

for dirname in list_dirs_data:
    data_samplename = dirname.split("crab_HAA4bAnalysis_")[1]
    print "Processing data sample dir " + dirname
    data_scale_factor = 1.0
    write_data_string = data_samplename + " " + str(data_scale_factor) + "\n"
    out_file.write(write_data_string)
out_file.close()

print "Data and Sample Normalizations saved in output file"
print "All done!"
