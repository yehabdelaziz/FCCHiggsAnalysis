import ROOT
import os
import argparse
import config_final
import json


def load_process_info(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data


def apply_cuts_and_fill_histograms(files):


    inputDir=config_final.input_dir
    outputDir=config_final.output_dir 

    if not os.path.exists(outputDir):
        os.makedirs(outputDir) 

    process_info = load_process_info("/cvmfs/fcc.cern.ch/FCCDicts/FCCee_procDict_winter2023_IDEA.json")

    for filename in files:
        cross_section = process_info[filename]['crossSection']
        Total_number_of_events = files[filename]
        print("*************************************************")
        print("Running on process:",filename)

       # Create RDataFrame from the TTree
        rdf = ROOT.RDataFrame("events", f"{inputDir}{filename}.root") 
        num_entries_before = rdf.Count().GetValue()
        print("Number of raw events before cuts:",num_entries_before)
        print("Number of normalized events before cuts:",num_entries_before*cross_section*config_final.Luminosity/Total_number_of_events)
        rdf = rdf.Filter(config_final.cut_condition)
        num_entries_after = rdf.Count().GetValue()
        print("Number of raw events after cuts:",num_entries_after)
        print("Number of normalized events after cuts:",num_entries_after*cross_section*config_final.Luminosity/Total_number_of_events)

        # Define the cut and histogram parameters
        histograms=[]
        for key,value in config_final.histoList.items():
            hist = rdf.Histo1D((key, "", value["bin"], value["xmin"],value["xmax"]), value["name"])
            hist.Scale(cross_section*config_final.Luminosity/Total_number_of_events)
            histograms.append(hist)
         
    # Optionally, save the histograms to ROOT files
        output_file = ROOT.TFile(f"{outputDir}{filename}_sel_hist.root", "RECREATE")

    # Define List of histograms and add them
        for histo in histograms:
            histo.Write()
        output_file.Close()



if __name__ == "__main__":
    file_list = config_final.files

    apply_cuts_and_fill_histograms(file_list)
