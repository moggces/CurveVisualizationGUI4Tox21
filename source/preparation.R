
###### load mapping file

#mapping_file <- './data/cas_tox21id_mapping.txt'
mapping_file <- './data/tox21_id_12807_cas_structure_4v4_080813_v4.txt'
mapping <- load_mapping_file(mapping_file)
# out_name <- 'tr_inh_main'
# cebs_file <- 'U:/Projects/TOX21/CEBS/tr_inhibition_curvep-wauc.txt'
# cebs <- load_cebs_file(cebs_file, mapping, pathway="tr_antagonism", readout="luc")
# save(cebs, file=paste("./data/", out_name, ".RData", sep=""))
# 
# out_name <- 'tr_inh_via'
# cebs_file <- 'U:/Projects/TOX21/CEBS/tr_viability_curvep-wauc.txt'
# cebs <- load_cebs_file(cebs_file, mapping, pathway="tr_antagonism", readout="via")
# save(cebs, file=paste("./data/", out_name, ".RData", sep=""))

# out_name <- 'p53_act_main'
# cebs_file <- 'U:/Projects/TOX21/CEBS/p53_ratio_activation_curvep-wauc.txt'
# cebs <- load_cebs_file(cebs_file, mapping, pathway="p53", readout="ratio")
# save(cebs, file=paste("./data/", out_name, ".RData", sep=""))


out_names <- c('atad5_act_main', 'atad5_act_via', 
               'p53_act_main', 'p53_act_via', 'p53_act_ch2', 'p53_act_ch1', 
               'srf_inh_main', 'dsb_inh_main', 'srf_inh_via', 'dsb_inh_via',
               'tr_inh_main', 'tr_inh_via',
               'fxr_act_main', 'fxr_act_via', 'fxr_act_ch2', 'fxr_act_ch1',
               'fxr_inh_main', 'fxr_inh_via', 'fxr_inh_ch2', 'fxr_inh_ch1'
               )
cebs_files <- c('atad5_agonist_curvep-wauc.txt', 'atad5_viability_curvep-wauc.txt', 
                'p53_ratio_activation_curvep-wauc.txt', 'p53_viability_curvep-wauc.txt', 'p53_ch2_activation_curvep-wauc.txt', 'p53_ch1_inhibition_curvep-wauc.txt', 
                'dt40_657_viability_curvep-wauc.txt', 'dt40_100_viability_curvep-wauc.txt', 'dt40_653_viability_curvep-wauc.txt','dt40_653_viability_curvep-wauc.txt',
                'tr_inhibition_curvep-wauc.txt', 'tr_viability_curvep-wauc.txt',
                'fxr_ratio_activation_curvep-wauc.txt', 'fxr_viability_activation_curvep-wauc.txt', 'fxr_ch2_activation_curvep-wauc.txt', 'fxr_ch1_activation_curvep-wauc.txt',
                'fxr_ratio_inhibition_curvep-wauc.txt', 'fxr_viability_inhibition_curvep-wauc.txt', 'fxr_ch2_inhibition_curvep-wauc.txt', 'fxr_ch1_inhibition_curvep-wauc.txt'
                )
pathways <- c('atad5', 'atad5', 'p53', 'p53', 'p53', 'p53', 
              'dna_damage_srf', 'dna_damage_dsb', 'dna_damage_srf', 'dna_damage_dsb',
              'tr_antagonism', 'tr_antagonism',
              'fxr_agonism','fxr_agonism','fxr_agonism','fxr_agonism',
              'fxr_antagonism', 'fxr_antagonism','fxr_antagonism','fxr_antagonism'
              )
readouts <- c('luc', 'via', 'ratio', 'via', 'ch2', 'ch1', 
              'luc', 'luc', 'via', 'via',
              'luc', 'via',
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1')

for (i in 1:length(out_names))
{
  out_name <- out_names[i]
  cebs_file <- paste("U:/Projects/TOX21/CEBS/", cebs_files[i], sep="")
  cebs <- load_cebs_file(cebs_file, mapping, pathway=pathways[i], readout=readouts[i])
  save(cebs, file=paste("./data/", out_name, ".RData", sep=""))
}


