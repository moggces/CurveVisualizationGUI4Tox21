
###### load mapping file

#mapping_file <- './data/cas_tox21id_mapping.txt'
mapping_file <- './data/tox21_id_12807_cas_structure_4v4_080813_v4.txt'
mapping <- load_mapping_file(mapping_file)

out_name <- 'autofluo_hek293_cell_blue'
cebs_file <- 'U:/Projects/TOX21/auto_fluro/hek293_autofluo_cell_blue_curvep_2.5sd.txt'
cebs <- load_cebs_file(cebs_file, mapping, pathway=NULL, readout=NULL)
save(cebs, file=paste("./data/", out_name, ".RData", sep=""))

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


# out_names <- c('atad5_act_main', 'atad5_act_via', 
#                'p53_act_main', 'p53_act_via', 'p53_act_ch2', 'p53_act_ch1', 
#                'srf_inh_main', 'dsb_inh_main', 'srf_inh_via', 'dsb_inh_via',
#                'tr_inh_main', 'tr_inh_via',
#                'fxr_act_main', 'fxr_act_via', 'fxr_act_ch2', 'fxr_act_ch1',
#                'fxr_inh_main', 'fxr_inh_via', 'fxr_inh_ch2', 'fxr_inh_ch1'
#                )

out_name <- 'tr_act_main'
cebs_file <- 'U:/Projects/TOX21/CEBS/tr_activation_curvep-wauc.txt'
cebs <- load_cebs_file(cebs_file, mapping, pathway="tr_agonism", readout="luc")
save(cebs, file=paste("./data/", out_name, ".RData", sep=""))

# cebs_files <- c('atad5_agonist_curvep-wauc.txt', 'atad5_viability_curvep-wauc.txt', 
#                 'p53_ratio_activation_curvep-wauc.txt', 'p53_viability_curvep-wauc.txt', 'p53_ch2_activation_curvep-wauc.txt', 'p53_ch1_inhibition_curvep-wauc.txt', 
#                 'dt40_657_viability_curvep-wauc.txt', 'dt40_100_viability_curvep-wauc.txt', 'dt40_653_viability_curvep-wauc.txt','dt40_653_viability_curvep-wauc.txt',
#                 'tr_inhibition_curvep-wauc.txt', 'tr_viability_curvep-wauc.txt',
#                 'fxr_ratio_activation_curvep-wauc.txt', 'fxr_viability_activation_curvep-wauc.txt', 'fxr_ch2_activation_curvep-wauc.txt', 'fxr_ch1_activation_curvep-wauc.txt',
#                 'fxr_ratio_inhibition_curvep-wauc.txt', 'fxr_viability_inhibition_curvep-wauc.txt', 'fxr_ch2_inhibition_curvep-wauc.txt', 'fxr_ch1_inhibition_curvep-wauc.txt'
#                 )
# pathways <- c('atad5', 'atad5', 'p53', 'p53', 'p53', 'p53', 
#               'dna_damage_srf', 'dna_damage_dsb', 'dna_damage_srf', 'dna_damage_dsb',
#               'tr_antagonism', 'tr_antagonism',
#               'fxr_agonism','fxr_agonism','fxr_agonism','fxr_agonism',
#               'fxr_antagonism', 'fxr_antagonism','fxr_antagonism','fxr_antagonism'
#               )

# readouts <- c('luc', 'via', 'ratio', 'via', 'ch2', 'ch1', 
#               'luc', 'luc', 'via', 'via',
#               'luc', 'via',
#               'ratio', 'via', 'ch2', 'ch1',
#               'ratio', 'via', 'ch2', 'ch1')

# 
# out_names <- c('pparg_act_main', 'pparg_act_ch2', 'pparg_act_ch1',
#                'are_act_main', 'are_act_via', 'are_act_ch2', 'are_act_ch1', 
#                'aromatase_inh_main', 'aromatase_inh_via'
# 
# cebs_files <- c('pparg_ratio_activation_curvep-wauc.txt', 'pparg_ch2_activation_curvep-wauc.txt', 'pparg_ch1_activation_curvep-wauc.txt',
#                 'are_ratio_activation_curvep-wauc.txt', 'are_viability_activation_curvep-wauc.txt', 'are_ch2_activation_curvep-wauc.txt', 'are_ch1_activation_curvep-wauc.txt',
#                 'aromatase_inhibition_curvep-wauc.txt', 'aromatase_viability_inhibition_curvep-wauc.txt'
#                 )
# 
# pathways <- c(rep('pparg_agonism', 3), 
#               rep('nrf2/are',4),
#               rep('aromatase_antagonism', 2)
#               )
# 
# 
# readouts <- c('ratio', 'ch2', 'ch1',
#               'ratio', 'via', 'ch2', 'ch1',
#               'luc', 'via'
#               )


out_names <- c('pparg_inh_main', 'pparg_inh_via','pparg_inh_ch2', 'pparg_inh_ch1',
               'ppard_act_main', 'ppard_act_via','ppard_act_ch2', 'ppard_act_ch1',
               'ppard_inh_main', 'ppard_inh_via','ppard_inh_ch2', 'ppard_inh_ch1',
               'gr_act_main', 'gr_act_ch2', 'gr_act_ch1',
               'gr_inh_main', 'gr_inh_via', 'gr_inh_ch2', 'gr_inh_ch1',
               'arfull_act_main', 'arfull_inh_main', 'arfull_inh_via',
               'erfull_act_main', 'erfull_inh_main', 'erfull_inh_via',
               'arpartial_act_main', 'arpartial_act_ch2', 'arpartial_act_ch1',
               'arpartial_inh_main','arpartial_inh_via', 'arpartial_inh_ch2', 'arpartial_inh_ch1',
               'erpartial_act_main', 'erpartial_act_ch2', 'erpartial_act_ch1',
               'erpartial_inh_main', 'erpartial_inh_via', 'erpartial_inh_ch2', 'erpartial_inh_ch1',
               'hse_act_main', 'hse_act_via','hse_act_ch2', 'hse_act_ch1',
               'mito_inh_main', 'mito_inh_via','mito_inh_ch2', 'mito_inh_ch1',
               'ahr_act_main', 'ahr_act_via'
)

cebs_files <- c('pparg_ratio_inhibition_curvep-wauc.txt', 'pparg_viability_inhibition_curvep-wauc.txt', 'pparg_ch2_inhibition_curvep-wauc.txt', 'pparg_ch1_inhibition_curvep-wauc.txt',
                'ppard_ratio_activation_curvep-wauc.txt', 'ppard_viability_activation_curvep-wauc.txt', 'ppard_ch2_activation_curvep-wauc.txt', 'ppard_ch1_activation_curvep-wauc.txt',
                'ppard_ratio_inhibition_curvep-wauc.txt', 'ppard_viability_inhibition_curvep-wauc.txt', 'ppard_ch2_inhibition_curvep-wauc.txt', 'ppard_ch1_inhibition_curvep-wauc.txt',
                'gr_ratio_activation_curvep-wauc.txt', 'gr_ch2_activation_curvep-wauc.txt', 'gr_ch1_activation_curvep-wauc.txt',
                'gr_ratio_inhibition_curvep-wauc.txt', 'gr_viability_inhibition_curvep-wauc.txt', 'gr_ch2_inhibition_curvep-wauc.txt', 'gr_ch1_inhibition_curvep-wauc.txt',
                'mdakb2ar_activation_curvep-wauc.txt', 'mdakb2ar_inhibition_curvep-wauc.txt', 'mdakb2ar_viability_curvep-wauc.txt',
                'bg1er_activation_curvep-wauc.txt', 'bg1er_inhibition_curvep-wauc.txt', 'bg1er_viability_curvep-wauc.txt',
                'hek293ar_ratio_activation_curvep-wauc.txt', 'hek293ar_ch2_activation_curvep-wauc.txt', 'hek293ar_ch1_activation_curvep-wauc.txt',
                'hek293ar_ratio_inhibition_curvep-wauc.txt', 'hek293ar_viability_inhibition_curvep-wauc.txt', 'hek293ar_ch2_inhibition_curvep-wauc.txt', 'hek293ar_ch1_inhibition_curvep-wauc.txt',
                'hek293er_ratio_activation_curvep-wauc.txt', 'hek293er_ch2_activation_curvep-wauc.txt', 'hek293er_ch1_activation_curvep-wauc.txt',
                'hek293er_ratio_inhibition_curvep-wauc.txt', 'hek293er_viability_inhibition_curvep-wauc.txt', 'hek293er_ch2_inhibition_curvep-wauc.txt', 'hek293er_ch1_inhibition_curvep-wauc.txt',
                'hse_ratio_activation_curvep-wauc.txt', 'hse_viability_activation_curvep-wauc.txt', 'hse_ch2_activation_curvep-wauc.txt', 'hse_ch1_activation_curvep-wauc.txt',
                'mito_ratio_inhibition_curvep-wauc.txt', 'mito_viability_inhibition_curvep-wauc.txt', 'mito_rho_inhibition_curvep-wauc.txt', 'mito_fitc_inhibition_curvep-wauc.txt',
                'ahr_activation_curvep-wauc.txt','ahr_viability_curvep-wauc.txt'
                )


pathways <- c(rep('pparg_antagonism', 4), 
              rep('ppard_agonism',4),
              rep('ppard_antagonism', 4), 
              rep('gr_agonism',3),
              rep('gr_antagonism', 4), 
              rep('ar_agonism(f)', 1), rep('ar_antagonism(f)', 2),
              rep('er_agonism(f)', 1), rep('er_antagonism(f)', 2),
              rep('ar_agonism(p)',3),
              rep('ar_antagonism(p)', 4),
              rep('er_agonism(p)',3),
              rep('er_antagonism(p)', 4),
              rep('hse', 4),
              rep('mitotox', 4),
              rep('ahr_agonism', 2)
)


readouts <- c(
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'luc', 'luc', 'via',
              'luc', 'luc', 'via',
              'ratio', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'ratio', 'via', 'ch2', 'ch1',
              'luc', 'via'
)
              

for (i in 1:length(out_names))
{
  out_name <- out_names[i]
  cebs_file <- paste("U:/Projects/TOX21/CEBS/", cebs_files[i], sep="")
  cebs <- load_cebs_file(cebs_file, mapping, pathway=pathways[i], readout=readouts[i])
  save(cebs, file=paste("./data/", out_name, ".RData", sep=""))
}


