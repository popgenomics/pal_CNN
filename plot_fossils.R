library(tidyverse)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)
simul_name = args[1]
work_dir = args[2]
setwd(work_dir)

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

infile = sprintf('simulated_fossils_%s.txt', simul_name)
pbdb = read.table(infile, h=T)
pbdb = as_tibble(pbdb)

pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)

pbdb = pbdb %>% group_by(genus, species) %>% summarise(simulation, geography ,max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()

pbdb = pbdb %>% dplyr::mutate(median = (max_ma+min_ma)/2)

pbdb = pbdb %>% dplyr::arrange(median)

pbdb = pbdb %>% dplyr::mutate(y=1:nrow(pbdb))

plot_count = pbdb %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y)) + geom_segment(size=3) + 
  xlim(0, max_time) + ylim(0, max_id) + theme_void()

ggsave(sprintf('%s_count.png',simul_name), width = 20, height = 20, units = "cm", bg='white')

# demes
infile = sprintf('simulated_fossils_%s.txt',simul_name)
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
  
  file_name = sprintf('%s_D%s.png', simul_name, deme)
  ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')
}
infile = sprintf('df_to_plot_%s.csv', simul_name)
data = read.csv(infile, h=T, sep = '\t')

plot_mspd = data %>% ggplot(aes(x=age_Ma, y=mspd, xend=age_Ma, yend=mspd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_sp_by_deme) + theme_void()

file_name = sprintf('%s_mspd.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')

plot_vspd = data %>% ggplot(aes(x=age_Ma, y=vspd, 
                                xend=age_Ma, yend=vspd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_sp_by_deme) + theme_void()

file_name = sprintf('%s_vspd.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')

plot_mgd = data %>% ggplot(aes(x=age_Ma, y=mgd, 
                               xend=age_Ma, yend=mgd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_genera_by_deme) + theme_void()

file_name = sprintf('%s_mgd.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')

plot_vgd = data %>% ggplot(aes(x=age_Ma, y=vgd, 
                               xend=age_Ma, yend=vgd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_genera_by_deme) + theme_void()

file_name = sprintf('%s_vgd.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')

plot_msp = data %>% ggplot(aes(x=age_Ma, y=msp, 
                               xend=age_Ma, yend=msp)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_sp_by_genera) + theme_void()

file_name = sprintf('%s_msp.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')

plot_vsp = data %>% ggplot(aes(x=age_Ma, y=vsp, 
                               xend=age_Ma, yend=vsp)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_sp_by_genera) + theme_void()

file_name = sprintf('%s_vsp.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')

plot_ng = data %>% ggplot(aes(x=age_Ma, y=ng, 
                              xend=age_Ma, yend=ng)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_genera) + theme_void()

file_name = sprintf('%s_ng.png',simul_name)
ggsave(file_name, width = 20, height = 20, units = "cm", bg='white')