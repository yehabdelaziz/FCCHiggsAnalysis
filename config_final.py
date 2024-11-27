input_dir="stage1/"
output_dir="final/"


cut_condition=" fourMuons_pmin > 5 && pmiss < 20 && vis_e_other_particles > 95   && non_res_Z_m < 65   && non_res_Z_m > 10"

Luminosity = 10800000 # in pb

histoList = {
"hist_fourMuons_mass":{"name":"fourMuons_mass","label":"M_{4#mu} [GeV]","bin":20,"xmin":120.0,"xmax":130.0},
"hist_Z_res_mass":{"name":"zll_mass","label":"Z mass [GeV]","bin":50,"xmin":0.0,"xmax":250.0},
"hist_Z_nonres_mass":{"name":"non_res_Z_m","label":"Z* mass [GeV]","bin":50,"xmin":0.0,"xmax":250.0},
"hist_vis_e_wo_muons":{"name":"vis_e_other_particles","label":"Visible Energy w/o muons [GeV]","bin":50,"xmin":0.0,"xmax":250.0},
"hist_missing_p":{"name":"pmiss","label":"missing energy","bin":50,"xmin":0.0,"xmax":250.0},

}



files = { 
'wzp6_ee_qqH_HZZ_llll_ecm240':1200000,
'p8_ee_ZZ_ecm240':56162093,
'p8_ee_WW_ecm240':373375386,
'p8_ee_Zqq_ecm240':100559248,


}
