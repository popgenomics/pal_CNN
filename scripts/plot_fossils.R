library(tidyverse)
# events
infile = 'simulated_mass_extinctions_S5.txt'
events = read.table(infile, h=T, sep='\t')
events = as_tibble(events)

# fossiles
infile = 'simulated_fossils_S5.txt'
pbdb = read.table(infile, h=T)
pbdb = as_tibble(pbdb)

#pbdb = pbdb %>% dplyr::sample_n(size=700000, replace=F)
pbdb = pbdb %>% dplyr::mutate(genus=as.factor(genus), species=as.factor(species), L=max_ma-min_ma)

pbdb = pbdb %>% group_by(genus, species) %>% summarise(simulation, max_ma=max(max_ma), min_ma=min(min_ma)) %>% ungroup()

pbdb = pbdb %>% dplyr::mutate(median = (max_ma+min_ma)/2)
#pbdb = pbdb %>% dplyr::arrange(genus)
pbdb = pbdb %>% dplyr::arrange(median)

pbdb = pbdb %>% dplyr::mutate(y=1:nrow(pbdb))
#x %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y, label=accepted_name)) + geom_segment(size=3) + geom_text(hjust='outside', nudge_x = -60) +

plot_1 = pbdb %>% ggplot(aes(x=max_ma, y=y, xend=min_ma, yend=y)) + geom_segment(size=3)  +
	scale_color_viridis_d() +
	xlab("Time (Ma)") + ylab('Species ID')

for(time in events$time){
	plot_1 = plot_1 + geom_vline(xintercept=time)
}

plot_1

