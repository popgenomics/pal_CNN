import pandas as pd
from numpy import var, mean

def setup_range(path, nb_simulations, scenario, fileout):
    """
    Read the specified number simulation in the input directory.
    Designed to work with simulation named "S{number}".

    Parameters
    ----------
    path : str
        The path of the directory with the simulations.
    nb_simulations : int
        The number of simulations to analyzes.
    fileout : str
        The path of the output directory where the file will be store.

    Returns
    -------
    None
    """
    time_max = 0
    nb_id_max = 0
    nb_id_max_deme = 0
    nb_genus_max = 0
    mean_genus_max_deme = 0
    var_genus_max_deme = 0
    mean_sp_max_deme = 0
    var_sp_max_deme = 0
    mean_sp_max_genera = 0
    var_sp_max_genera = 0
    
    for i in range(nb_simulations):
        file = f'{path}/df_to_plot_S{scenario}_{i}.csv'
        df = pd.read_csv(file, sep='\t', header=0)
        for time in range(600+1):
            current_line = df.iloc[[time]]
            time_max = max(time_max, current_line['age_Ma'].iloc[0])
            nb_id_max = max(nb_id_max, current_line['nb_id_max'].iloc[0])
            nb_id_max_deme = max(nb_id_max_deme, current_line['nb_id_max_deme'].iloc[0])
            nb_genus_max = max(nb_genus_max, current_line['ng'].iloc[0])
            mean_genus_max_deme = max(mean_genus_max_deme, current_line['mgd'].iloc[0])
            var_genus_max_deme = max(var_genus_max_deme, current_line['vgd'].iloc[0])
            mean_sp_max_deme = max(mean_sp_max_deme, current_line['mspd'].iloc[0])
            var_sp_max_deme = max(var_sp_max_deme, current_line['vspd'].iloc[0])
            mean_sp_max_genera = max(mean_sp_max_genera, current_line['msp'].iloc[0])
            var_sp_max_genera = max(var_sp_max_genera, current_line['vsp'].iloc[0])
    
    with open(fileout, 'w') as fileOut:
        fileOut.write('time_max\tnb_id_max\tnb_id_max_deme\tmean_sp_max_genera\tvar_sp_max_genera\tmean_sp_max_deme\tvar_sp_max_deme\tmean_genus_max_deme\tvar_genus_max_deme\tnb_genus_max\n')
        fileOut.write(f'{time_max}\t{nb_id_max}\t{nb_id_max_deme}\t{mean_sp_max_genera}\t{var_sp_max_genera}\t{mean_sp_max_deme}\t{var_sp_max_deme}\t{mean_genus_max_deme}\t{var_genus_max_deme}\t{nb_genus_max}\n')

