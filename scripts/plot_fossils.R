library(tidyverse)
library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly = TRUE)
simul_name = args[1]
txt_dir = args[2]
csv_dir = args[3]
range_graph = args[4]
size_graph = args[5]
graph_dir = args[6]
nb_demes = args[7]
nb_demes = as.integer(nb_demes)
size_graph = as.integer(size_graph)
setwd(txt_dir)

parameters = read.table(range_graph, h=T, sep = '\t')
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
max_sdt = parameters['standard_deviation_time'][1,1]
max_sdtd = parameters['standard_deviation_time_deme'][1,1]
max_number_sp = parameters['number_species'][1,1]
max_number_spd = parameters['number_species_deme'][1,1]
nb_demes = nb_demes -1

infile = sprintf('simulated_fossils_%s.txt', simul_name)
pbdb = read.table(infile, h=T)
pbdb = as_tibble(pbdb)

pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)

pbdb = pbdb %>% group_by(genus, species) %>% summarise(simulation ,max_ma=max(max_ma), min_ma=min(min_ma)) %>% dplyr::distinct() %>% ungroup()

pbdb = pbdb %>% dplyr::mutate(median = (max_ma+min_ma)/2)

pbdb = pbdb %>% dplyr::arrange(median)

pbdb = pbdb %>% dplyr::mutate(y=1:nrow(pbdb))

plot_median = pbdb %>% ggplot(aes(x=median, y=y)) + xlim(0, max_time) + ylim(0, max_id) + geom_line()+ theme_void()

ggsave(sprintf('%s_median.png',simul_name), width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_count = pbdb %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y)) + geom_segment(size=3) + 
  xlim(0, max_time) + ylim(0, max_id) + theme_void()

ggsave(sprintf('%s_count.png',simul_name), width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

# demes
infile = sprintf('simulated_fossils_%s.txt',simul_name)
pbdb = read.table(infile, h=T)
pbdb = as_tibble(pbdb)
pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)
demes = pbdb %>% select('geography') %>% distinct()
demes = demes %>% dplyr::arrange(geography)
for(deme in 0:nb_demes){
  pbdb_deme = pbdb %>% filter(geography == deme)
  pbdb_deme = pbdb_deme %>% group_by(genus, species) %>% summarise(simulation, geography ,max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()
  pbdb_deme = pbdb_deme %>% dplyr::mutate(median = (max_ma+min_ma)/2)
  pbdb_deme = pbdb_deme %>% dplyr::arrange(median) %>% distinct()
  pbdb_deme = pbdb_deme %>% dplyr::mutate(y=row_number())
  name = sprintf('plot_deme_%s', deme)
  plot_deme = pbdb_deme %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y)) + geom_segment(size=3)  +
    scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_id_by_deme) + theme_void()
  assign(name, plot_deme)
  
  file_name = sprintf('%s_D%s.png', simul_name, deme)
  ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)
  
  plot_deme_median = pbdb_deme %>% ggplot(aes(x=median, y=y)) + xlim(0, max_time) + ylim(0, max_id_by_deme) + geom_line()+ theme_void()
  file_name = sprintf('%s_D%s_median.png', simul_name, deme)
  ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)
}

setwd(csv_dir)

infile = sprintf('%s.csv', simul_name)
data = read.csv(infile, h=T, sep = '\t')

data = data %>% replace(is.na(.), 0)

plot_mspd = data %>% ggplot(aes(x=age_Ma, y=mspd, xend=age_Ma, yend=mspd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_sp_by_deme) + theme_void()

file_name = sprintf('%s_mspd.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_vspd = data %>% ggplot(aes(x=age_Ma, y=vspd, 
                                xend=age_Ma, yend=vspd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_sp_by_deme) + theme_void()

file_name = sprintf('%s_vspd.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_mgd = data %>% ggplot(aes(x=age_Ma, y=mgd, 
                               xend=age_Ma, yend=mgd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_genera_by_deme) + theme_void()

file_name = sprintf('%s_mgd.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_vgd = data %>% ggplot(aes(x=age_Ma, y=vgd, 
                               xend=age_Ma, yend=vgd)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_genera_by_deme) + theme_void()

file_name = sprintf('%s_vgd.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_msp = data %>% ggplot(aes(x=age_Ma, y=msp, 
                               xend=age_Ma, yend=msp)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_mean_sp_by_genera) + theme_void()

file_name = sprintf('%s_msp.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_vsp = data %>% ggplot(aes(x=age_Ma, y=vsp, 
                               xend=age_Ma, yend=vsp)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_var_sp_by_genera) + theme_void()

file_name = sprintf('%s_vsp.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_ng = data %>% ggplot(aes(x=age_Ma, y=ng, 
                              xend=age_Ma, yend=ng)) + geom_line()  +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_genera) + theme_void()

file_name = sprintf('%s_ng.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdt = data %>% ggplot(aes(x=age_Ma, y = sdt, xend=age_Ma, yend=sdt)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdt) + theme_void()

file_name = sprintf('%s_sdt.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd0 = data %>% ggplot(aes(x=age_Ma, y = sdtd0, xend=age_Ma, yend=sdtd0)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd0.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd1 = data %>% ggplot(aes(x=age_Ma, y = sdtd1, xend=age_Ma, yend=sdtd1)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd1.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd2 = data %>% ggplot(aes(x=age_Ma, y = sdtd2, xend=age_Ma, yend=sdtd2)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd2.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd3 = data %>% ggplot(aes(x=age_Ma, y = sdtd3, xend=age_Ma, yend=sdtd3)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd3.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd4 = data %>% ggplot(aes(x=age_Ma, y = sdtd4, xend=age_Ma, yend=sdtd4)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd4.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd5 = data %>% ggplot(aes(x=age_Ma, y = sdtd5, xend=age_Ma, yend=sdtd5)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd5.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd6 = data %>% ggplot(aes(x=age_Ma, y = sdtd6, xend=age_Ma, yend=sdtd6)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd6.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd7 = data %>% ggplot(aes(x=age_Ma, y = sdtd7, xend=age_Ma, yend=sdtd7)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd7.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd8 = data %>% ggplot(aes(x=age_Ma, y = sdtd8, xend=age_Ma, yend=sdtd8)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd8.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_sdtd9 = data %>% ggplot(aes(x=age_Ma, y = sdtd9, xend=age_Ma, yend=sdtd9)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_sdtd) + theme_void()

file_name = sprintf('%s_sdtd9.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nsp = data %>% ggplot(aes(x=age_Ma, y = nsp, xend=age_Ma, yend=nsp)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_sp) + theme_void()

file_name = sprintf('%s_nsp.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd0 = data %>% ggplot(aes(x=age_Ma, y = nspd0, xend=age_Ma, yend=nspd0)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd0.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd1 = data %>% ggplot(aes(x=age_Ma, y = nspd1, xend=age_Ma, yend=nspd1)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd1.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd2 = data %>% ggplot(aes(x=age_Ma, y = nspd2, xend=age_Ma, yend=nspd2)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd2.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd3 = data %>% ggplot(aes(x=age_Ma, y = nspd3, xend=age_Ma, yend=nspd3)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd3.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd4 = data %>% ggplot(aes(x=age_Ma, y = nspd4, xend=age_Ma, yend=nspd4)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd4.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd5 = data %>% ggplot(aes(x=age_Ma, y = nspd5, xend=age_Ma, yend=nspd5)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd5.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd6 = data %>% ggplot(aes(x=age_Ma, y = nspd6, xend=age_Ma, yend=nspd6)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd6.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd7 = data %>% ggplot(aes(x=age_Ma, y = nspd7, xend=age_Ma, yend=nspd7)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd7.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd8 = data %>% ggplot(aes(x=age_Ma, y = nspd8, xend=age_Ma, yend=nspd8)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd8.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)

plot_nspd9 = data %>% ggplot(aes(x=age_Ma, y = nspd9, xend=age_Ma, yend=nspd9)) + geom_line() +
  scale_color_viridis_d() + xlim(0, max_time) + ylim(0, max_number_spd) + theme_void()

file_name = sprintf('%s_nspd9.png',simul_name)
ggsave(file_name, width = size_graph, height = size_graph, units = "cm", bg='white', path = graph_dir)