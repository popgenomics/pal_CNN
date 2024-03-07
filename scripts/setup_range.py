import pandas as pd
from numpy import var, mean

def setup_range(path, nb_simulations, fileout):
    time_max = None
    nb_id_max = 0
    nb_id_max_deme = 0
    nb_genus_max = 0
    mean_genus_max_deme = 0
    var_genus_max_deme = 0
    mean_sp_max_deme = 0
    var_sp_max_deme = 0
    mean_sp_max_genera = 0
    var_sp_max_genera = 0
    
    for i in range(1,nb_simulations+1):
        file = f'{path}/simulated_fossils_S{i}.txt'
        array = pd.read_csv(file, sep='\t', header=0)
        # Count max number of id and max number of genera
        if i == 1 :
            time_max = array['max_ma'].max(axis=0)
            nb_id_max = array.shape[0]
        else:
            time_max = max(time_max, array['max_ma'].max(axis=0))
            nb_id_max = max(nb_id_max, array.shape[0])
            
        # Count max number of id in deme
        for j in range(array['geography'].nunique()):
            deme_array = array.loc[array['geography']==j]             
            if nb_id_max_deme <= deme_array.shape[0]:
                nb_id_max_deme = deme_array.shape[0]
        
        for time in range(time_max+1):
            array_at_time = array[(array['max_ma'] >= time) & (array['min_ma'] <= time)]
            nb_genus_max = max(nb_genus_max, array_at_time['genus'].nunique())
            nb_sp_deme = []
            nb_genera_deme = []
            # max mean and variance number of species / genera within demes over time
            for j in array_at_time['geography'].unique():
                array_deme_time =  array_at_time[array_at_time['geography']==j]
                nb_sp_deme.append(array_deme_time['species'].count())
                nb_genera_deme.append(array_deme_time['genus'].nunique())

            if mean_sp_max_deme <= mean(nb_sp_deme):
                mean_sp_max_deme = mean(nb_sp_deme)
            if len(nb_sp_deme) > 1 and var_sp_max_deme <= var(nb_sp_deme, ddof=1):
                var_sp_max_deme = var(nb_sp_deme, ddof=1)
            
            if mean_genus_max_deme <= mean(nb_genera_deme):
                mean_genus_max_deme = mean(nb_genera_deme)
            if len(nb_genera_deme) > 1  and var_genus_max_deme <= var(nb_genera_deme, ddof=1):
                var_genus_max_deme = var(nb_genera_deme, ddof=1)
            
            nb_sp_genus = []
            array_at_time_simple = array_at_time.drop(['simulation', 'geography', 'specimen', 'extinction_statu', 'max_ma', 'min_ma'], axis=1)
            # max mean and variance number of species within genera over time
            for j in array_at_time_simple['genus'].unique():
                array_genus_time = array_at_time_simple[array_at_time_simple['genus']==j]
                nb_sp_genus.append(array_genus_time['species'].nunique())
            if mean_sp_max_genera <= mean(nb_sp_genus):
                mean_sp_max_genera = mean(nb_sp_genus)
            if len(nb_sp_genus) > 1 and var_sp_max_genera <= var(nb_sp_genus, ddof=1):
                var_sp_max_genera = var(nb_sp_genus, ddof=1)
    
    with open(fileout, 'w') as fileOut:
        fileOut.write('time_max\tnb_id_max\tnb_id_max_deme\tmean_sp_max_genera\tvar_sp_max_genera\tmean_sp_max_deme\tvar_sp_max_deme\tmean_genus_max_deme\tvar_genus_max_deme\tnb_genus_max\n')
        fileOut.write(f'{time_max}\t{nb_id_max}\t{nb_id_max_deme}\t{mean_sp_max_genera}\t{var_sp_max_genera}\t{mean_sp_max_deme}\t{var_sp_max_deme}\t{mean_genus_max_deme}\t{var_genus_max_deme}\t{nb_genus_max}\n')

