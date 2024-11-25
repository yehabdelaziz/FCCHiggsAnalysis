#ifndef ZHfunctions_H
#define ZHfunctions_H

#include <cmath>
#include <vector>
#include <math.h>

#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "ReconstructedParticle2MC.h"
//#include <algorithm> 
namespace FCCAnalyses { namespace ZHfunctions {

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> getTwoHighestPMuons(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons) {
    if (muons.size() < 2) {
        std::cerr << "Not enough muons to select the top 2." << std::endl;
        return {}; // Return an empty vector if there are fewer than 2 muonsy
    }
    // Create a copy of the input vector to sort
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> sortedMuons = muons;

    // Lambda function to sort muons by pt in descending order
    auto compareByPt = [&](edm4hep::ReconstructedParticleData a, edm4hep::ReconstructedParticleData b) {
        return sqrt(a.momentum.x * a.momentum.x + a.momentum.y*a.momentum.y + a.momentum.z*a.momentum.z  ) > sqrt(b.momentum.x * b.momentum.x + b.momentum.y*b.momentum.y + b.momentum.z*b.momentum.z );
    };

    // Sort the copied vector by pt in descending order
    std::sort(sortedMuons.begin(), sortedMuons.end(), compareByPt);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selected_muons;
 
    for (const auto& muon : sortedMuons) {
        if (selected_muons.empty()) {
            selected_muons.push_back(muon);
        } else if (selected_muons[0].charge != muon.charge) {
            selected_muons.push_back(muon);
            break;
        }
    }
    
    if (selected_muons.size() == 2) {
        return selected_muons;
    } 
    else {
	return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>(); // Return an empty vector if no valid pair is found
       }
    // Return the result vector
}




ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> getTwoHighestPtMuons(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> muons) {
    if (muons.size() < 2) {
        std::cerr << "Not enough muons to select the top 2." << std::endl;
        return {}; // Return an empty vector if there are fewer than 2 muons
    }
    // Create a copy of the input vector to sort
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> sortedMuons = muons;

    // Lambda function to sort muons by pt in descending order
    auto compareByPt = [&](edm4hep::ReconstructedParticleData a, edm4hep::ReconstructedParticleData b) {
        return sqrt(a.momentum.x * a.momentum.x + a.momentum.y*a.momentum.y) > sqrt(b.momentum.x * b.momentum.x + b.momentum.y*b.momentum.y);
    };

    // Sort the copied vector by pt in descending order
    std::sort(sortedMuons.begin(), sortedMuons.end(), compareByPt);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> selected_muons;
    for (const auto& muon : sortedMuons) {
        if (selected_muons.empty()) {
            selected_muons.push_back(muon);
        } else if (selected_muons[0].charge != muon.charge) {
            selected_muons.push_back(muon);
            break;
        }
    }

    if (selected_muons.size() == 2) {
        return selected_muons;
    }
    else {
        return ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>(); // Return an empty vector if no valid pair is found
       }
    // Return the result vector
}


ROOT::VecOps::RVec<TLorentzVector>  makeLorentzVectors(ROOT::VecOps::RVec<float>  jets_px, ROOT::VecOps::RVec<float>  jets_py, ROOT::VecOps::RVec<float>  jets_pz, ROOT::VecOps::RVec<float>  jets_e){
ROOT::VecOps::RVec<TLorentzVector> result;
    for(int i=0; i<jets_px.size(); i++) {
        TLorentzVector tlv;
        tlv.SetPxPyPzE(jets_px[i], jets_py[i], jets_pz[i], jets_e[i]);
        result.push_back(tlv);
    }
    return result;


   }

float InvariantMass(const TLorentzVector &tlv1, const TLorentzVector &tlv2)
    {
      float E = tlv1.E() + tlv2.E();
      float px = tlv1.Px() + tlv2.Px();
      float py = tlv1.Py() + tlv2.Py();
      float pz = tlv1.Pz() + tlv2.Pz();
      return std::sqrt(E * E - px * px - py * py - pz * pz);
    }



// compute the cone isolation for reco particles
struct coneIsolation0 {

    coneIsolation0(float arg_dr_min, float arg_dr_max);
    double deltaR(double eta1, double phi1, double eta2, double phi2) { return TMath::Sqrt(TMath::Power(eta1-eta2, 2) + (TMath::Power(phi1-phi2, 2))); };

    float dr_min = 0;
    float dr_max = 0.4;
    ROOT::VecOps::RVec<float>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rps) ;
};

coneIsolation0::coneIsolation0(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
ROOT::VecOps::RVec<float>  coneIsolation0::coneIsolation0::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rps) {

    ROOT::VecOps::RVec<float>  result;
    result.reserve(in.size());

    std::vector<ROOT::Math::PxPyPzEVector> lv_reco;
    std::vector<ROOT::Math::PxPyPzEVector> lv_charged;
    std::vector<ROOT::Math::PxPyPzEVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);
        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }

    for(size_t i = 0; i < in.size(); ++i) {

        ROOT::Math::PxPyPzEVector tlv;
        tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
        lv_reco.push_back(tlv);
    }


    // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
    for (auto & lv_reco_ : lv_reco) { 
        double sumNeutral = 0.0;
        double sumCharged = 0.0;

        // charged
        for (auto & lv_charged_ : lv_charged) {

            double dr = coneIsolation0::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_charged_.Eta(), lv_charged_.Phi());
            if(dr >= dr_min && dr < dr_max) sumCharged += lv_charged_.P();
        }

        // neutral
        for (auto & lv_neutral_ : lv_neutral) {

            double dr = coneIsolation0::deltaR(lv_reco_.Eta(), lv_reco_.Phi(), lv_neutral_.Eta(), lv_neutral_.Phi());
            if(dr >= dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
        }

        double sum = sumCharged + sumNeutral;
        double ratio= sum / lv_reco_.P();
        result.emplace_back(ratio);
    }
    return result;
}

// build the Z resonance based on the available leptons. Returns the best lepton pair compatible with the Z mass and recoil at 125 GeV
// technically, it returns a ReconstructedParticleData object with index 0 the di-lepton system, index and 2 the leptons of the pair
struct resonanceBuilder_mass {
    float m_resonance_mass;
    bool m_use_MC_Kinematics;
    resonanceBuilder_mass(float arg_resonance_mass, bool arg_use_MC_Kinematics);
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents, ROOT::VecOps::RVec<int> daugthers) ;
};

resonanceBuilder_mass::resonanceBuilder_mass(float arg_resonance_mass, bool arg_use_MC_Kinematics) {m_resonance_mass = arg_resonance_mass, m_use_MC_Kinematics = arg_use_MC_Kinematics;}

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> resonanceBuilder_mass::resonanceBuilder_mass::operator()(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> legs, ROOT::VecOps::RVec<int> recind, ROOT::VecOps::RVec<int> mcind, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> reco, ROOT::VecOps::RVec<edm4hep::MCParticleData> mc, ROOT::VecOps::RVec<int> parents, ROOT::VecOps::RVec<int> daugthers) {

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    result.reserve(3);
    std::vector<std::vector<int>> pairs; // for each permutation, add the indices of the muons
    int n = legs.size();
 
    if(n > 1) {
        ROOT::VecOps::RVec<bool> v(n);
        std::fill(v.end() - 2, v.end(), true); // helper variable for permutations
        do {
            std::vector<int> pair;
            edm4hep::ReconstructedParticleData reso;
            reso.charge = 0;
            TLorentzVector reso_lv; 
            for(int i = 0; i < n; ++i) {
                if(v[i]) {
                    pair.push_back(i);
                    reso.charge += legs[i].charge;
                    TLorentzVector leg_lv;

                    if(m_use_MC_Kinematics) { // MC kinematics
                        int track_index = legs[i].tracks_begin;   // index in the Track array
                        int mc_index = ReconstructedParticle2MC::getTrack2MC_index(track_index, recind, mcind, reco);
                        if (mc_index >= 0 && mc_index < mc.size()) {
                            leg_lv.SetXYZM(mc.at(mc_index).momentum.x, mc.at(mc_index).momentum.y, mc.at(mc_index).momentum.z, mc.at(mc_index).mass);
                        }
                    }
                    else { // reco kinematics
                         leg_lv.SetXYZM(legs[i].momentum.x, legs[i].momentum.y, legs[i].momentum.z, legs[i].mass);
                    }

                    reso_lv += leg_lv;
                }
            }

            if(reso.charge != 0) continue; // neglect non-zero charge pairs
            reso.momentum.x = reso_lv.Px();
            reso.momentum.y = reso_lv.Py();
            reso.momentum.z = reso_lv.Pz();
            reso.mass = reso_lv.M();
            result.emplace_back(reso);
            pairs.push_back(pair);

        } while(std::next_permutation(v.begin(), v.end()));
    }
    else {
        std::cout << "ERROR: resonanceBuilder_mass, at least two leptons required." << std::endl;
        //exit(1);
    }

    if(result.size() > 1) {
        ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> bestReso;

        int idx_min = -1;
        float d_min = 9e9;
        for (int i = 0; i < result.size(); ++i) {

           // float d = std::pow(result.at(i).mass - m_resonance_mass, 2); // mass
            float d = abs(result.at(i).mass - m_resonance_mass); // mass
            if(d < d_min) {
                d_min = d;
                idx_min = i;
            }
        }
        if(idx_min > -1) { 
            bestReso.push_back(result.at(idx_min));
            auto & l1 = legs[pairs[idx_min][0]];
            auto & l2 = legs[pairs[idx_min][1]];
            bestReso.emplace_back(l1);
            bestReso.emplace_back(l2);
        }
        else {
            std::cout << "ERROR: resonanceBuilder_mass, no mininum found." << std::endl;
            exit(1);
        }
        return bestReso;
    }
    else {
        auto & l1 = legs[0];
        auto & l2 = legs[1];
        result.emplace_back(l1);
        result.emplace_back(l2);
        return result;
    }
}




struct sel_iso {
    sel_iso(float arg_max_iso);
    float m_max_iso = .25;
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<float>  iso);
  };

sel_iso::sel_iso(float arg_max_iso) : m_max_iso(arg_max_iso) {};
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>  sel_iso::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<float>  iso) {
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> result;
    result.reserve(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        auto & p = in[i];
        if (iso[i] < m_max_iso) {
            result.emplace_back(p);
        }
    }
    return result;
}

 
// compute the cone isolation for reco particles
struct coneIsolation {

    coneIsolation(float arg_dr_min, float arg_dr_max);
    float dr_min = 0;
    float dr_max = 0.4;
    ROOT::VecOps::RVec<float>  operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rps) ;
};

coneIsolation::coneIsolation(float arg_dr_min, float arg_dr_max) : dr_min(arg_dr_min), dr_max( arg_dr_max ) { };
ROOT::VecOps::RVec<float>  coneIsolation::coneIsolation::operator() (ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> rps) {

    ROOT::VecOps::RVec<float>  result;
    result.reserve(in.size());

    ROOT::VecOps::RVec<TLorentzVector> lv_reco;
    ROOT::VecOps::RVec<TLorentzVector> lv_charged;
    ROOT::VecOps::RVec<TLorentzVector> lv_neutral;

    for(size_t i = 0; i < rps.size(); ++i) {
        TLorentzVector tlv;
        tlv.SetPxPyPzE(rps.at(i).momentum.x, rps.at(i).momentum.y, rps.at(i).momentum.z, rps.at(i).energy);

        if(rps.at(i).charge == 0) lv_neutral.push_back(tlv);
        else lv_charged.push_back(tlv);
    }

    for(size_t i = 0; i < in.size(); ++i) {
        TLorentzVector tlv;
        tlv.SetPxPyPzE(in.at(i).momentum.x, in.at(i).momentum.y, in.at(i).momentum.z, in.at(i).energy);
        lv_reco.push_back(tlv);
    }

    // compute the isolation (see https://github.com/delphes/delphes/blob/master/modules/Isolation.cc#L154) 
    for (auto & lv_reco_ : lv_reco) {
        double sumNeutral = 0.0;
        double sumCharged = 0.0;
        // charged
        for (auto & lv_charged_ : lv_charged) {
            double dr = lv_reco_.Vect().Angle(lv_charged_.Vect());
            if(dr > dr_min && dr < dr_max) sumCharged += lv_charged_.P();
        }

        // neutral
        for (auto & lv_neutral_ : lv_neutral) {
            double dr = lv_reco_.Vect().Angle(lv_neutral_.Vect());
            if(dr > dr_min && dr < dr_max) sumNeutral += lv_neutral_.P();
        }
        double sum = sumCharged + sumNeutral;
        double ratio= sum / lv_reco_.P();
        result.emplace_back(ratio);
    }
    return result;
}
 
 
 
// returns missing energy vector, based on reco particles
ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> missingEnergy(float ecm, ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> in, float p_cutoff = 0.0) {
    float px = 0, py = 0, pz = 0, e = 0;
    for(auto &p : in) {
        if (std::sqrt(p.momentum.x * p.momentum.x + p.momentum.y*p.momentum.y) < p_cutoff) continue;
        px += -p.momentum.x;
        py += -p.momentum.y;
        pz += -p.momentum.z;
        e += p.energy;
    }

    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> ret;
    edm4hep::ReconstructedParticleData res;
    res.momentum.x = px;
    res.momentum.y = py;
    res.momentum.z = pz;
    res.energy = ecm-e;
    ret.emplace_back(res);
    return ret;
}

// calculate the cosine(theta) of the missing energy vector
float get_cosTheta_miss(ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> met){
    float costheta = 0.;
    if(met.size() > 0) {
        TLorentzVector lv_met;
        lv_met.SetPxPyPzE(met[0].momentum.x, met[0].momentum.y, met[0].momentum.z, met[0].energy);
        costheta = fabs(std::cos(lv_met.Theta()));
    }
    return costheta;
}



}}

#endif
