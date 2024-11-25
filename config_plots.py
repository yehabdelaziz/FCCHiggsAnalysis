import ROOT

histograms_to_plot = ["hist_fourMuons_mass","hist_Z_res_mass","hist_Z_nonres_mass","hist_vis_e_wo_muons","hist_missing_p"]
input_file_directory="./final/"
output_dir="./plots"


background_prefix = {
'WW':'p8_ee_WW_ecm240',
'Zqq':'p8_ee_Zqq_ecm240',
'ZZ':'p8_ee_ZZ_ecm240',
}

signal_prefix = {
'qqH_HZZ':'wzp6_ee_qqH_HZZ_llll_ecm240',
}

legends = {}
legends['ZZ'] = 'ZZ'
legends['Zqq'] = 'Zqq'
legends['WW'] = 'WW'

legends['qqH_HZZ']='Z(qq)H(4#mu)'

colors={}
colors['ZZ'] = ROOT.kBlue
colors['Zqq'] = ROOT.kRed
colors['WW'] = ROOT.kGray

colors['qqH_HZZ']= ROOT.kPink









