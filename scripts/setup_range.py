import pandas as pd
import os

def setup_range(path, fileout, nb_time):
    """
    Read the specified number simulation in the input directory.

    Parameters
    ----------
    path : str
        The path of the directory with the simulations.
    fileout : str
        The path of the output directory where the file will be store.
    nb_time : int
        The number of time in simulation.
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
    standard_deviation_time = 0
    standard_deviation_time_deme = 0
    # standard_deviation_time_deme0 = 0
    # standard_deviation_time_deme1 = 0
    # standard_deviation_time_deme2 = 0
    # standard_deviation_time_deme3 = 0
    # standard_deviation_time_deme4 = 0
    # standard_deviation_time_deme5 = 0
    # standard_deviation_time_deme6 = 0
    # standard_deviation_time_deme7 = 0
    # standard_deviation_time_deme8 = 0
    # standard_deviation_time_deme9 = 0
    # standard_deviation_time_deme = 0
    number_species = 0
    number_species_deme = 0
    # number_species_deme0 = 0
    # number_species_deme1 = 0
    # number_species_deme2 = 0
    # number_species_deme3 = 0
    # number_species_deme4 = 0
    # number_species_deme5 = 0
    # number_species_deme6 = 0
    # number_species_deme7 = 0
    # number_species_deme8 = 0
    # number_species_deme9 = 0


    files_list = [path + '/' + f for f in os.listdir(path)]
    for file in files_list:  # for each files
        df = pd.read_csv(file, sep='\t', header=0)
        for time in range(nb_time+1): # for each times
            current_line = df.iloc[[time]]
            # Check the maximum value for each stats
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
            standard_deviation_time = max(standard_deviation_time, current_line['sdt'].iloc[0])
            standard_deviation_time_deme = max(standard_deviation_time_deme, current_line['sdtd0'].iloc[0],
                                               current_line['sdtd1'].iloc[0], current_line['sdtd2'].iloc[0],
                                               current_line['sdtd3'].iloc[0], current_line['sdtd4'].iloc[0],
                                               current_line['sdtd5'].iloc[0], current_line['sdtd6'].iloc[0],
                                               current_line['sdtd7'].iloc[0], current_line['sdtd8'].iloc[0],
                                               current_line['sdtd9'].iloc[0])
            # standard_deviation_time_deme0 = max(standard_deviation_time_deme0, current_line['sdtd0'].iloc[0])
            # standard_deviation_time_deme1 = max(standard_deviation_time_deme1, current_line['sdtd1'].iloc[0])
            # standard_deviation_time_deme2 = max(standard_deviation_time_deme2, current_line['sdtd2'].iloc[0])
            # standard_deviation_time_deme3 = max(standard_deviation_time_deme3, current_line['sdtd3'].iloc[0])
            # standard_deviation_time_deme4 = max(standard_deviation_time_deme4, current_line['sdtd4'].iloc[0])
            # standard_deviation_time_deme5 = max(standard_deviation_time_deme5, current_line['sdtd5'].iloc[0])
            # standard_deviation_time_deme6 = max(standard_deviation_time_deme6, current_line['sdtd6'].iloc[0])
            # standard_deviation_time_deme7 = max(standard_deviation_time_deme7, current_line['sdtd7'].iloc[0])
            # standard_deviation_time_deme8 = max(standard_deviation_time_deme8, current_line['sdtd8'].iloc[0])
            # standard_deviation_time_deme9 = max(standard_deviation_time_deme9, current_line['sdtd9'].iloc[0])
            number_species = max(number_species, current_line['nsp'].iloc[0])
            number_species_deme = max(number_species_deme, current_line['nspd0'].iloc[0], current_line['nspd1'].iloc[0],
                                      current_line['nspd2'].iloc[0], current_line['nspd3'].iloc[0], current_line['nspd4'].iloc[0],
                                      current_line['nspd5'].iloc[0], current_line['nspd6'].iloc[0], current_line['nspd7'].iloc[0],
                                      current_line['nspd8'].iloc[0], current_line['nspd9'].iloc[0])
            # number_species_deme0 = max(number_species_deme0, current_line['nspd0'].iloc[0])
            # number_species_deme1 = max(number_species_deme1, current_line['nspd1'].iloc[0])
            # number_species_deme2 = max(number_species_deme2, current_line['nspd2'].iloc[0])
            # number_species_deme3 = max(number_species_deme3, current_line['nspd3'].iloc[0])
            # number_species_deme4 = max(number_species_deme4, current_line['nspd4'].iloc[0])
            # number_species_deme5 = max(number_species_deme5, current_line['nspd5'].iloc[0])
            # number_species_deme6 = max(number_species_deme6, current_line['nspd6'].iloc[0])
            # number_species_deme7 = max(number_species_deme7, current_line['nspd7'].iloc[0])
            # number_species_deme8 = max(number_species_deme8, current_line['nspd8'].iloc[0])
            # number_species_deme9 = max(number_species_deme9, current_line['nspd9'].iloc[0])

    # create the output
    with open(fileout, 'w') as fileOut:
        fileOut.write('time_max\tnb_id_max\tnb_id_max_deme\tmean_sp_max_genera\tvar_sp_max_genera\tmean_sp_max_deme'
                      '\tvar_sp_max_deme\tmean_genus_max_deme\tvar_genus_max_deme\tnb_genus_max\tstandard_deviation_time'
                      '\tstandard_deviation_time_deme\tnumber_species\tnumber_species_deme\n')
        fileOut.write(f'{time_max}\t{nb_id_max}\t{nb_id_max_deme}\t{mean_sp_max_genera}\t{var_sp_max_genera}\t'
                      f'{mean_sp_max_deme}\t{var_sp_max_deme}\t{mean_genus_max_deme}\t{var_genus_max_deme}\t{nb_genus_max}'
                      f'\t{standard_deviation_time}\t{standard_deviation_time_deme}\t{number_species}'
                      f'\t{number_species_deme}\n')

