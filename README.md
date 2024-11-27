# Higgs to 4 leptons analysis at the future circular collider
 In this analysis we study e+ e- -> ZH events where Z-> jj and H -> ZZ* -> 4 leptons. The samples used for this analysis are available 
in [FCC-ee Delphes winter 2023 samples](https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/).
This analysis is divided into three parts: stage1 - final - plot.

In stage1, events are read from the Delphes root file and the 4 leptons are chosen from the list of Reconstructed particle objects.
An event selection is applied to choose events containing at least 4 leptons.
The Higgs is reconstructed from the 4 leptons as follows:
First pair of leptons (From On-shell Z)​:
- Oppositely charged leptons
- The pair which minimises |M<sub>ll</sub> - M<sub>Z</sub>|​.<br>

Second Pair of leptons (From off-shell Z)​:
- Oppositely charged leptons​
- Highest momentum pair of the remaining leptons​.<br>

In final stage, analysis cuts are applied on the ROOT tree produced from stage1 and histgorams are created from the produced branches. 
The histograms are then weighted according to each process's cross section and total integrated luminosity.
The plot stage is where the signal and backgrounds are plotted. The background is stacked and the signal is overlayed above the background.
To run the analysis, first source the [key4hep](https://github.com/key4hep/) stack.
```
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2024-03-10
```
To run stage1, 

```
python3 analysis_stage1.py process_name number_of_events
``` 

where process_name is the name of the process name found [here](https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/)
and the number_of_events is the number of events you wish to run.<br>

To run final step, first open the file config_final.py specify the cut conditions,luminosity, histograms and the list of processes.
Then run:

```
python3 analysis_final.py
```

To run the plot stage, open the config_plots file and specify the list of histograms you want to plot, the list of signal and background
processes and the plot legend and colors for each process. Then run:


```
python3 plotting_FCC.py
```

