import ROOT 
import os
import argparse



ROOT.gInterpreter.Declare('#include "src/ReconstructedParticle.cc"')
ROOT.gInterpreter.Declare('#include "src/functions.h"')

parser = argparse.ArgumentParser(description="Higgs Analysis at FCC")

parser.add_argument("process", help="Name of the process")
parser.add_argument("n_events", help="Number of events you want to process",type=int)

args = parser.parse_args()


output_dir = "stage1"
current_directory = os.getcwd()
new_directory_path = os.path.join(current_directory, output_dir)

if not os.path.exists(new_directory_path):
    os.makedirs(new_directory_path)

# Get absolute paths of all files
file_directory = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/{}".format(args.process)
file_list = [os.path.abspath(os.path.join(file_directory, f)) for f in os.listdir(file_directory) if os.path.isfile(os.path.join(file_directory, f))]
print(args.n_events)


def analyse(df):

    df = df.Alias("Particle0", "Particle#0.index")
    df = df.Alias("Particle1", "Particle#1.index")
    df = df.Alias("MCRecoAssociations0", "MCRecoAssociations#0.index")
    df = df.Alias("MCRecoAssociations1", "MCRecoAssociations#1.index")


    df = df.Alias('Muon0', 'Muon#0.index')
    df = df.Define('muons','ReconstructedParticle::get(Muon0, ReconstructedParticles)')
    df = df.Define('muons_p','ReconstructedParticle::get_p(muons)')

    df = df.Define('selected_muons','ReconstructedParticle::sel_p(2)(muons)')
    df = df.Define('selected_muons_n','ReconstructedParticle::get_n(selected_muons)')
    df = df.Define('selected_muons_p','ReconstructedParticle::get_p(selected_muons)')
    df = df.Filter('selected_muons_n > 3')

            # Find muons pair with mass closest to 91.2 (On-shell Z)
    df = df.Define("zbuilder_result", "FCCAnalyses::ZHfunctions::resonanceBuilder_mass(91.2,false)(selected_muons, MCRecoAssociations0, MCRecoAssociations1, ReconstructedParticles, Particle, Particle0, Particle1)")

    df = df.Define("zll", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> {zbuilder_result[0]}") 
    df = df.Define("zll_charge", 'ReconstructedParticle::get_charge(zll)') 
    df = df.Define("zll_mass", "ReconstructedParticle::get_mass(zll)[0]") 


    df = df.Define("zll_muons", "ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> {zbuilder_result[1],zbuilder_result[2]}") # the leptons
    df = df.Define('zll_muons_p','ReconstructedParticle::get_p(zll_muons)')
    df = df.Define('zll_muons_n','ReconstructedParticle::get_n(zll_muons)')
    df = df.Define('rest_of_muons',"ReconstructedParticle::remove(selected_muons,zll_muons)") #Remove the muon pair of the on-shell Z from the muon >


    df = df.Define("non_res_Z","FCCAnalyses::ZHfunctions::getTwoHighestPMuons(rest_of_muons)") # Find the higest p muon pair from the remaining muon>
    df = df.Define('non_res_Z_p','ReconstructedParticle::get_p(non_res_Z)')
    df = df.Define('non_res_Z_n','ReconstructedParticle::get_n(non_res_Z)')
    df = df.Define('non_res_Z_px','ReconstructedParticle::get_px(non_res_Z)')
    df = df.Define('non_res_Z_py','ReconstructedParticle::get_py(non_res_Z)')
    df = df.Define('non_res_Z_pz','ReconstructedParticle::get_pz(non_res_Z)')
    df = df.Define('non_res_Z_e','ReconstructedParticle::get_e(non_res_Z)')
    df = df.Define('non_res_Z_tlv','FCCAnalyses::ZHfunctions::makeLorentzVectors(non_res_Z_px,non_res_Z_py,non_res_Z_pz,non_res_Z_e)')
    df = df.Define('non_res_Z_m','FCCAnalyses::ZHfunctions::InvariantMass(non_res_Z_tlv[0],non_res_Z_tlv[1])')      
    df = df.Define('non_res_Z_angle','non_res_Z_tlv[0].Vect().Angle(non_res_Z_tlv[1].Vect())')
 
    df = df.Define('fourMuons',"ReconstructedParticle::merge(zll_muons,non_res_Z)") #Merge the two muon pairs
    df = df.Define('fourMuons_p',"ReconstructedParticle::get_p(fourMuons)") 
    df = df.Define('fourMuons_pmin',"return *std::min_element(fourMuons_p.begin(), fourMuons_p.end());") 
    df = df.Define('fourMuons_p4','ReconstructedParticle::get_P4vis(fourMuons)') # P4 of four muon pairs
    df = df.Define('fourMuons_mass','fourMuons_p4.M()') # Mass of 4 muon pairs
    df = df.Define("rest_of_particles","ReconstructedParticle::remove(ReconstructedParticles,fourMuons)")

    df = df.Define("vis_p4_other_particles","ReconstructedParticle::get_P4vis(rest_of_particles)") 
    df = df.Define("vis_e_other_particles","vis_p4_other_particles.E()") 


    df = df.Define("Emiss","FCCAnalyses::ZHfunctions::missingEnergy(240,ReconstructedParticles)") 
    df = df.Define("pmiss","Emiss[0].energy")
    df = df.Define("cosTheta_miss", "FCCAnalyses::ZHfunctions::get_cosTheta_miss(Emiss)") 
    df = df.Define("fourMuons_iso","FCCAnalyses::ZHfunctions::coneIsolation(0.0,0.523599)(fourMuons,rest_of_particles)") 
    df = df.Define('fourMuons_min_iso',"return *std::max_element(fourMuons_iso.begin(),fourMuons_iso.end());") 
    return df



branchlist = [
                      "selected_muons_n",
                      "selected_muons_p",
                      "fourMuons_p",
                      "fourMuons_mass",
                      "zll_mass",
                      "non_res_Z_m",
                      "vis_e_other_particles",
                      "fourMuons_pmin",
                      "non_res_Z_angle",
                      "fourMuons_iso",
                      "fourMuons_min_iso",
                      "pmiss",
                      "cosTheta_miss",
                      ]

if __name__ == "__main__":

    df = ROOT.RDataFrame("events", file_list)
    df = df.Range(args.n_events)
    df = analyse(df)
    df.Snapshot("events", "stage1/{}.root".format(args.process), branchlist)
