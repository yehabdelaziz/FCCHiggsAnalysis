# Higgs to 4 leptons analysis at the future circular collider
 In this analysis we study e+ e- -> ZH events where Z-> jj and H -> ZZ* -> 4 leptons. The samples used for this analysis are available 
in [FCC-ee Delphes winter 2023 samples](https://fcc-physics-events.web.cern.ch/fcc-ee/delphes/winter2023/idea/)
This analysis is divided into three parts: stage1 - final - plot. In stage1, events are read from the Delphes root file and the 4 leptons are
chosen from the list of Reconstructed particle objects. An event selection is applied to choose events containing at least 4 leptons.
The Higgs is reconstructed from the 4 leptons as follows:
First pair of leptons (From On-shell Z)​:
- Oppositely charged leptons
- The pair which minimises |M<sub>ll</sub> - M<sub>Z</sub>|​.
Second Pair of leptons (From off-shell Z)​:
- Oppositely charged leptons​
- Highest momentum pair of the remaining leptons​. 
