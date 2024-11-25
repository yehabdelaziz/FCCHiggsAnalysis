import ROOT
from ROOT import TFile, THStack, TCanvas, TLegend, gStyle,TLatex
import json
import math
import numpy as np
import config_final
import config_plots
# Dictionaries for background and signal processes
ROOT.gROOT.SetBatch(True)

legends = config_plots.legends
colors = config_plots.colors

input_file_directory = config_plots.input_file_directory

def load_m_vis_histogram(file, legend_text, hist_name):
    hist = None
    key = file.Get(hist_name)
    if key:
        hist = key.Clone()
        hist.SetDirectory(0)  # Detach histogram from file
        hist.SetTitle(legend_text)
    return hist




# Main function to plot histograms
def plot_histograms(hist_name):
    # Create a canvas
    canvas = TCanvas("canvas", "Histograms", 800, 600)

    # Create a stack for background histograms
    stack = THStack("stack", "FCCAnalyses: FCC-ee Simulation (Delphes)")

    legend = TLegend(0.6, 0.7, 0.9, 0.9)  
        
    # List to store signal histograms
    signal_histograms = []
    background_histograms = []
    sum_bckg=0
    # Load background histograms
    print("Loading background processes...")
    for process, prefix in config_plots.background_prefix.items():
        file_path = f"{input_file_directory}/{prefix}_sel_hist.root"
        file = TFile.Open(file_path)
        if file and not file.IsZombie():
            hist = load_m_vis_histogram(file, legends[process],hist_name)
            if hist:
                hist.SetLineColor(colors[process])
                hist.SetFillColor(colors[process])
                background_histograms.append(hist)  # Append histogram to signal list
                print(prefix,":",hist.Integral())
                sum_bckg+=hist.Integral()
                stack.Add(hist)
                legend.AddEntry(hist, legends[process], "f")  # Add line for signal histograms
  # Add histogram to stack
            file.Close()
    print("Total Background:",sum_bckg)


    # Load signal histograms
    for process, prefix in config_plots.signal_prefix.items():
        file_path = f"{input_file_directory}/{prefix}_sel_hist.root"
        file = TFile.Open(file_path)
        if file and not file.IsZombie():
            hist = load_m_vis_histogram(file, legends[process],hist_name)
            if hist:
                hist.SetLineColor(colors[process])
                signal_histograms.append(hist)  # Append histogram to signal list
                legend.AddEntry(hist, legends[process], "l")  # Add line for signal histograms
                print("signal",prefix,":",hist.Integral())
            file.Close()


    canvas.SetLogy()
    # Draw stacked background histograms
    stack.Draw("HIST")
    legend.Draw()
    # Draw signal histograms overlayed
    for hist in signal_histograms:
        hist.Draw("HIST SAME")



    # Save canvas as PDF
        # Adjust y-axis range
#    stack.SetMinimum(1e-1)
#    stack.SetMaximum(100)
    stack.GetXaxis().SetTitle(config_final.histoList[hist_name]["label"])
#    stack.GetXaxis().SetRangeUser(120,130)

    stack.GetYaxis().SetTitle("Events / 1 GeV")
    canvas.SaveAs(f"{config_plots.output_dir}/{histogram_name}.png")


if __name__ == "__main__":
    
    for histogram_name in config_plots.histograms_to_plot:
        plot_histograms(histogram_name)
