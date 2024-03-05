library(tidyverse)
library(ggplot2)

setwd("//wsl.localhost/Ubuntu/home/arthur_boddaert/pal_CNN/results")

infile = 'range_graph.txt'
parameters = read.table(infile, h=T, sep = '\t')
max_time = parameters['time_max'][1,1]
max_id = parameters['nb_id_max'][1,1]
max_id_by_deme = parameters['nb_id_max_deme'][1,1]
max_mean_sp_by_genera = parameters['mean_sp_max_genera'][1,1]
max_var_sp_by_genera = parameters['var_sp_max_genera'][1,1]
max_mean_sp_by_deme = parameters['mean_sp_max_deme'][1,1]
max_var_sp_by_deme = parameters['var_sp_max_deme'][1,1]
max_mean_genera_by_deme = parameters['mean_genus_max_deme'][1,1]
max_var_genera_by_deme = parameters['var_genus_max_deme'][1,1]
max_genera = parameters['nb_genus_max'][1,1]

for(i in c(1:3)){
  # events
  infile = sprintf('simulated_mass_extinctions_S%s.txt',i)
  events = read.table(infile, h=T, sep='\t')
  events = as_tibble(events)
  
  # fossiles
  infile = sprintf('simulated_fossils_S%s.txt',i)
  pbdb = read.table(infile, h=T)
  pbdb = as_tibble(pbdb)
  
  #pbdb = pbdb %>% dplyr::sample_n(size=700000, replace=F)
  pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)
  
  pbdb = pbdb %>% group_by(genus, species) %>% summarise(simulation, geography ,max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
  
  pbdb = pbdb %>% dplyr::mutate(median = (max_ma+min_ma)/2)
  #pbdb = pbdb %>% dplyr::arrange(genus)
  pbdb = pbdb %>% dplyr::arrange(median)
  
  pbdb = pbdb %>% dplyr::mutate(y=1:nrow(pbdb))
  #x %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y, label=accepted_name)) + geom_segment(size=3) + geom_text(hjust='outside', nudge_x = -60) +
  
  plot_count = pbdb %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y)) + geom_segment(size=3) + 
    xlim(0, max_time) + ylim(0, max_id) + theme_void()
  
  #for(time in events$time){
  #	plot_1 = plot_1 + geom_vline(xintercept=time)
  #}
  ggsave(sprintf('S%i_count.png',i), width = 20, height = 20, units = "cm", bg='white')
  
  # demes
  infile = sprintf('simulated_fossils_S%s.txt',i)
  pbdb = read.table(infile, h=T)
  pbdb = as_tibble(pbdb)
  pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)
  demes = pbdb %>% select('geography') %>% distinct()
  demes = demes %>% dplyr::arrange(geography)
  for(deme in demes$geography){
    pbdb_deme = pbdb %>% filter(geography == deme)
    pbdb_deme = pbdb_deme %>% group_by(genus, species) %>% summarise(simulation, geography ,max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
    pbdb_deme = pbdb_deme %>% dplyr::mutate(median = (max_ma+min_ma)/2)
    pbdb_deme = pbdb_deme %>% dplyr::arrange(median) %>% distinct()
    pbdb_deme = pbdb_deme %>% dplyr::mutate(y=row_number())
    plot_deme = pbdb_deme %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y)) + geom_segment(size=3)  +
      scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_id_by_deme) + theme_void()
    
    file_name = sprintf('S%s_D%s.png',i, deme)
    ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  }
  
  # mean number of species within demes 
  infile = sprintf('simulated_fossils_S%s.txt',i)
  pbdb = read.table(infile, h=T)
  pbdb = as_tibble(pbdb)
  pbdb$geography = as.factor(pbdb$geography)
  pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)
  pbdb = pbdb %>% group_by(geography, genus, species) %>% summarise(simulation, max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
  pbdb = pbdb %>% dplyr::distinct()
  max_nb_species = dim(pbdb)[1]
  
  age_Ma = c()
  mean_nb_species = c()
  mean_nb_genus = c()
  var_nb_species = c()
  var_nb_genus = c()
  for(j in c(seq(600, 0, by=-1))){
    species_at_time = pbdb %>% filter(max_ma >= j & min_ma <= j)
    genus_at_time = pbdb %>% select(!(species)) %>% filter(max_ma >= j & min_ma <= j)
    genus_at_time = genus_at_time %>% group_by(geography, genus) %>% summarise(simulation, max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
    genus_at_time = genus_at_time %>% unique()
    species_by_demes = fct_count(species_at_time$geography)
    genus_by_demes = fct_count(genus_at_time$geography)
    age_Ma = c(age_Ma, j)
    mean_nb_species = c(mean_nb_species, mean(species_by_demes$n))
    var_nb_species = c(var_nb_species, var(species_by_demes$n))
    mean_nb_genus = c(mean_nb_genus, mean(genus_by_demes$n))
    var_nb_genus = c(var_nb_genus, var(genus_by_demes$n))
  }
  evol_species_by_deme = data.frame(age_Ma, mean_nb_species, var_nb_species, mean_nb_genus, var_nb_genus)
  
  plot_mspd = evol_species_by_deme %>% ggplot(aes(x=age_Ma, y=mean_nb_species, 
                                                  xend=age_Ma, yend=mean_nb_species)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_sp_by_deme) + theme_void()
  
  file_name = sprintf('S%s_mspd.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  
  plot_vspd = evol_species_by_deme %>% ggplot(aes(x=age_Ma, y=var_nb_species, 
                                                  xend=age_Ma, yend=max_var_sp_by_deme)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_sp_by_deme) + theme_void()
  
  file_name = sprintf('S%s_vspd.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  
  plot_mgd = evol_species_by_deme %>% ggplot(aes(x=age_Ma, y=mean_nb_genus, 
                                                 xend=age_Ma, yend=mean_nb_genus)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_genera_by_deme) + theme_void()
  
  file_name = sprintf('S%s_mgd.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  
  plot_vgd = evol_species_by_deme %>% ggplot(aes(x=age_Ma, y=var_nb_genus, 
                                                 xend=age_Ma, yend=var_nb_genus)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_genera_by_deme) + theme_void()
  
  file_name = sprintf('S%s_vgd.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  
  # mean number of species within genus 
  infile = sprintf('simulated_fossils_S%s.txt',i)
  pbdb = read.table(infile, h=T)
  pbdb = as_tibble(pbdb)
  pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)
  pbdb = pbdb %>% group_by(genus, species) %>% summarise(simulation, max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
  pbdb = pbdb %>% dplyr::distinct()
  
  age_Ma = c()
  mean_nb_species = c()
  var_nb_species = c()
  for(j in c(seq(600, 0, by=-1))){
    species_at_time = pbdb %>% filter(max_ma >= j & min_ma <= j)
    species_by_genus = fct_count(species_at_time$genus)
    nb_species = data.frame(species_by_genus) %>% filter(species_by_genus$n != 0)
    age_Ma = c(age_Ma, j)
    mean_nb_species = c(mean_nb_species, mean(nb_species$n))
    var_nb_species = c(var_nb_species, var(nb_species$n))
  }
  evol_species_by_genus = data.frame(age_Ma, mean_nb_species, var_nb_species)
  plot_msp = evol_species_by_genus %>% ggplot(aes(x=age_Ma, y=mean_nb_species, 
                                                  xend=age_Ma, yend=mean_nb_species)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_sp_by_genera) + theme_void()
  
  file_name = sprintf('S%s_msp.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  
  plot_vsp = evol_species_by_genus %>% ggplot(aes(x=age_Ma, y=var_nb_species, 
                                                 xend=age_Ma, yend=var_nb_species)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_sp_by_genera) + theme_void()
  
  file_name = sprintf('S%s_vsp.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
  
  # number of genera over time
  infile = sprintf('simulated_fossils_S%s.txt',i)
  pbdb = read.table(infile, h=T)
  pbdb = as_tibble(pbdb)
  pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)
  pbdb = pbdb %>% group_by(genus) %>% summarise(simulation, max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
  pbdb = pbdb %>% dplyr::distinct()
  age_Ma = c()
  nb_genus = c()
  for(j in c(seq(600, 0, by=-1))){
    nb_genus_at_time = dim(pbdb %>%filter(max_ma >= j & min_ma <= j))[1]
    age_Ma = c(age_Ma, j)
    nb_genus = c(nb_genus, nb_genus_at_time)
  }
  evol_genus_by_time = data.frame(age_Ma, nb_genus)
  plot_ng = evol_genus_by_time %>% ggplot(aes(x=age_Ma, y=nb_genus, 
                                              xend=age_Ma, yend=nb_genus)) + geom_line()  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_genera) + theme_void()
  
  file_name = sprintf('S%s_ng.png',i)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
}