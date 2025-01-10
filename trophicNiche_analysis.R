#'""" Code for 'Trophic niche dynamics of two fish mesoconsumers in adjacent 
#'    coastal habitats with varying nutrient regimes'
#'    @authors: Rolando O. Santos, W. Ryan James, Jennifer S. Rehage, 
#'    Cody W. Eggenberger, Justin S. Lesser, Christopher J. Madden
#'    Date: 1/10/25 """

# Mixing models ----
library(MixSIAR)
options(max.print = 6000000)

# ACS Individual nested model with 4 prey groups 
mix = load_mix_data(file("data/AC_SnookTarponIso.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c("Species", "ind"),
                    fac_random=c(F, T),
                    fac_nested=c(F, T),
                    cont_effects=NULL)


source = load_source_data(file("data/sCELA_AC4prey.csv"),
                          source_factors=NULL,
                          conc_dep=F,
                          data_type="means",
                          mix)

discr = load_discr_data(file("data/CELA_TEF.csv"), mix)

# Make an isospace plot
#plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.ac = run_model(run="long", mix, source, discr, model_filename,
                    alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_ac = list(summary_save = TRUE,
                 summary_name = "data/AC4prey_ss",
                 sup_post = FALSE,
                 plot_post_save_pdf = FALSE,
                 plot_post_name = "lower_posterior_density",
                 sup_pairs = FALSE,
                 plot_pairs_save_pdf = FALSE,
                 plot_pairs_name = "lower_pairs_plot",
                 sup_xy = TRUE,
                 plot_xy_save_pdf = FALSE,
                 plot_xy_name = "lower_xy_plot",
                 gelman = TRUE,
                 heidel = FALSE,
                 geweke = TRUE,
                 diag_save = TRUE,
                 diag_name = "data/AC4prey_diagnostics",
                 indiv_effect = FALSE,
                 plot_post_save_png = F,
                 plot_pairs_save_png = FALSE,
                 plot_xy_save_png = FALSE)

output_JAGS(jags.ac, mix, source, output_ac)



# MCS Individual nested model with 4 prey groups 
mix = load_mix_data(file("data/MC_SnookTarponIso.csv"),
                    iso_names=c("d13C","d15N","d34S"),
                    factors= c("Species", "ind"),
                    fac_random=c(F, T),
                    fac_nested=c(F, T),
                    cont_effects=NULL)


source = load_source_data(file("data/sCELA_MC4prey.csv"),
                          source_factors=NULL,
                          conc_dep=F,
                          data_type="means",
                          mix)

discr = load_discr_data(file("data/CELA_TEF.csv"), mix)

# Make an isospace plot
#plot_data(filename="isospace_plot", plot_save_pdf=FALSE, plot_save_png=FALSE, mix,source,discr)

# Write the JAGS model file
model_filename = "MixSIAR_model.txt"
resid_err = T
process_err = T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

jags.mc = run_model(run="long", mix, source, discr, model_filename,
                    alpha.prior = 1, resid_err, process_err)

# Process JAGS output
output_mc = list(summary_save = TRUE,
                 summary_name = "data/mc4prey_ss",
                 sup_post = FALSE,
                 plot_post_save_pdf = FALSE,
                 plot_post_name = "lower_posterior_density",
                 sup_pairs = FALSE,
                 plot_pairs_save_pdf = FALSE,
                 plot_pairs_name = "lower_pairs_plot",
                 sup_xy = TRUE,
                 plot_xy_save_pdf = FALSE,
                 plot_xy_name = "lower_xy_plot",
                 gelman = TRUE,
                 heidel = FALSE,
                 geweke = TRUE,
                 diag_save = TRUE,
                 diag_name = "data/mc4prey_diagnostics",
                 indiv_effect = FALSE,
                 plot_post_save_png = F,
                 plot_pairs_save_png = FALSE,
                 plot_xy_save_png = FALSE)

output_JAGS(jags.mc, mix, source, output_mc)


# hypervolume analysis ----

# function to bootstrap hypervolumes to generate confidence intervals
# data = dataframe or tibble of z-scored values to make hvs
# num = number of bootstrapped hvs to make
# frac = fraction of data points to use to make bootstrapped hvs
#    ** can also be input as a whole number for how many  
#       data points to use in resampling of bootstraps
# Returns vector of bootstrapped hvs
HVboot = function(data, num, frac, vol = NULL){
  require(tidyverse)
  df = tibble(Volume = rep(NA, len = num))
  n = ceiling(frac*nrow(data))
  if(frac>1){
    n = frac
  }
  for (i in 1:num){
    rand = sample_n(data, n, replace = TRUE)
    require(hypervolume)
    randh = hypervolume_gaussian(rand,
                                 samples.per.point = 5000,
                                 kde.bandwidth = estimate_bandwidth(rand), 
                                 sd.count = 3, 
                                 quantile.requested = 0.95, 
                                 quantile.requested.type = "probability", 
                                 chunk.size = 1000, 
                                 verbose = F)
    df$Volume[i] = as.numeric(get_volume(randh))
    if (i/num == .1){
      cat('10% done\n')
    } else if (i/num == 0.25){
      cat('25% done\n')
    }else if (i/num == 0.50){
      cat('50% done\n')
    }else if (i/num == 0.75){
      cat('75% done\n')
    }else if (i/num == 0.9){
      cat('90% done\n')
    }
  }
  if (is.null(vol) == F){
    cat('Volume (95% CI) = ', vol, ' (', 
        quantile(df$Volume, 0.025, na.rm = T), '-', 
        quantile(df$Volume, 0.975, na.rm = T), ')\n')
  }
  return(df)
}


# load z scored mixing model data
df = read_csv('data/SnookTarpon_zHV.csv') 

# Snook in ACS
acSN = df %>% filter(Species == 'Snook', System == 'ACS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
acSNhv = hypervolume_gaussian(acSN, name = 'Snook ACS',
                              #samples.per.point = ceiling((10^(3 + sqrt(ncol(acSN))))/nrow(acSN)),
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(acSN), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 121.38
get_volume(acSNhv)

# calculate CI
# Volume (95% CI) =  121.4  ( 64.50002 - 206.9598 )       
acSNhvCI = HVboot(acSN, 100, 0.75, vol = 121.4)
acSNhvCI$Species = 'Snook'
acSNhvCI$System = 'ACS'


# Snook in MCS
mcSN = df %>% filter(Species == 'Snook', System == 'MCS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
mcSNhv = hypervolume_gaussian(mcSN, name = 'Snook MCS',
                              #samples.per.point = ceiling((10^(3 + sqrt(ncol(mcSN))))/nrow(mcSN)),
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(mcSN), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 3.04
get_volume(mcSNhv)

# calculate CI
# Volume (95% CI) =  3  ( 0.809157 - 6.022489 )       
mcSNhvCI = HVboot(mcSN, 100, 0.75, vol = 3.0)
mcSNhvCI$Species = 'Snook'
mcSNhvCI$System = 'MCS'

# Tarpon in ACS
acTR = df %>% filter(Species == 'Tarpon', System == 'ACS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
acTRhv = hypervolume_gaussian(acTR, name = 'Tarpon AC',
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(acTR), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 32.1
get_volume(acTRhv)

# calculate CI
# Volume (95% CI) =  32.1  ( 4.668775 - 72.18363 )     
acTRhvCI = HVboot(acTR, 100, 0.75, vol = 32.1)
acTRhvCI$Species = 'Tarpon'
acTRhvCI$System = 'ACS'

# Tarpon in MCS
mcTR = df %>% filter(Species == 'Tarpon', System == 'MCS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
mcTRhv = hypervolume_gaussian(mcTR, name = 'Tarpon MCS',
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(mcTR), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 1.23
get_volume(mcTRhv)

# calculate CI
# Volume (95% CI) =  1.2  ( 0.0753357 - 3.86419 )     
mcTRhvCI = HVboot(mcTR, 100, 0.75, vol = 1.2)
mcTRhvCI$Species = 'Tarpon'
mcTRhvCI$System = 'MCS'

# bind all data together for plotting 
hvCI = bind_rows(acSNhvCI, mcSNhvCI, acTRhvCI, mcTRhvCI)
write_csv(hvCI, 'data/hv_CI.csv')

# overlap 
# data1 = z scored data for hv 1
# data2 = z scored data for hv 2
# num = number of iterations to run for bootstrapping
# frac = fraction of data points to use to make bootstrapped hvs
#    ** can also be input as a whole number for how many  
#       data points to use in resampling of bootstraps
# returns df with overlap and amount unique of each hv
bootOverlap = function(data1,data2,num,frac){
  require(tidyverse)
  df = tibble(sorenson = rep(NA, num), unique1 = NA, unique2 = NA)
  n1 = ceiling(frac*nrow(data1))
  if(frac>1){
    n1 = frac
  }
  n2 = ceiling(frac*nrow(data2))
  if(frac>1){
    n2 = frac
  }
  for (i in 1:num){
    rand1 = sample_n(data1, n1, replace = T)
    rand2 = sample_n(data2, n2, replace = T)
    require(hypervolume)
    randh1 = hypervolume_gaussian(rand1,
                                  samples.per.point = 5000,
                                  kde.bandwidth = estimate_bandwidth(rand1), 
                                  sd.count = 3, 
                                  quantile.requested = 0.95, 
                                  quantile.requested.type = "probability", 
                                  chunk.size = 1000, 
                                  verbose = F)
    randh2 = hypervolume_gaussian(rand2,
                                  samples.per.point = 5000,
                                  kde.bandwidth = estimate_bandwidth(rand2), 
                                  sd.count = 3, 
                                  quantile.requested = 0.95, 
                                  quantile.requested.type = "probability", 
                                  chunk.size = 1000, 
                                  verbose = F)
    set = hypervolume_set(randh1,randh2,check.memory = F, verbose = F)
    df$sorenson[i] = hypervolume_overlap_statistics(set)[2]
    df$unique1[i] = hypervolume_overlap_statistics(set)[3]
    df$unique2[i] = hypervolume_overlap_statistics(set)[4]
    
    if (i/num == .1){
      cat('10% done\n')
    } else if (i/num == 0.25){
      cat('25% done\n')
    }else if (i/num == 0.50){
      cat('50% done\n')
    }else if (i/num == 0.75){
      cat('75% done\n')
    }else if (i/num == 0.9){
      cat('90% done\n')
    }
    
  }
  return(df)
}

# system comparisons
# snook
hvSN= hypervolume_join(acSNhv, mcSNhv)
SN = hypervolume_set(acSNhv, mcSNhv,  check.memory = F)
hypervolume_overlap_statistics(SN)

ov1 = bootOverlap(acSN, mcSN, 100, 0.75)
ov1$comp = 'System'
ov1$type = 'Snook'

# Tarpon
hvTR= hypervolume_join(acTRhv, mcTRhv)
TR = hypervolume_set(acTRhv, mcTRhv, check.memory = F)

hypervolume_overlap_statistics(TR)

ov2 = bootOverlap(acTR, mcTR, 100, 0.75)
ov2$comp = 'System'
ov2$type = 'Tarpon'

# species comparisons
# Alligator
hvAC= hypervolume_join(acSNhv, acTRhv)
AC = hypervolume_set(acSNhv, acTRhv,  check.memory = F)
hypervolume_overlap_statistics(AC)

ov3 = bootOverlap(acSN, acTR, 100, 0.75)
ov3$comp = 'Species'
ov3$type = 'ACS'

# McCormick
hvMC= hypervolume_join(mcSNhv, mcTRhv)
MC = hypervolume_set(mcSNhv, mcTRhv, check.memory = F, verbose = F)

hypervolume_overlap_statistics(MC)

ov4 = bootOverlap(mcSN, mcTR, 100, 0.75)
ov4$comp = 'Species'
ov4$type = 'MCS'

# bind data together
OV = bind_rows(ov1, ov2, ov3, ov4)
write_csv(OV, 'data/Ov_CI.csv')


# overlap metrics

# hypervolume_overlap_statistics(SN)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.02065068    0.04046571    0.97925928    0.17370518

# hypervolume_overlap_statistics(TR)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.02106570    0.04126218    0.97857912    0.44043478 

# hypervolume_overlap_statistics(AC)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.07033968    0.13143431    0.91689751    0.68586981 

# hypervolume_overlap_statistics(MC)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.07389729    0.13762451    0.90342465    0.76063307 

