import os.path

import pandas as pd
from numpy import var, mean, nan, std
from numpy.random import choice

pd.options.mode.chained_assignment = None  # default='warn'

def data_to_dataframes(data, nb_time, simul, fileout, nb_deme):
    """
    Process the data to create file easily read by R to plot the informations.

    Parameters
    ----------
    data : str
        The path of the file produce by the simulator.
    nb_time : int
        The number of time in simulation.
    simul : int
        The simulation number
    fileout : int
        The path of tje output file produce by this function.
    nb_deme : int
        The number of deme in simulation.
    Returns
    -------
    None
    """
    # read the data
    df = pd.read_csv(data, sep='\t', header=0)

    # Number of unique species during simulation
    df_genera_sp = df[['genus','species']]
    nb_sp = df_genera_sp.groupby('genus').nunique().sum().values[0]
    nb_id_max = [nb_sp] + [0] * nb_time
    nb_id_max_deme = 0

    # Max number of unique species in a deme during simulation
    for j in range(df['geography'].nunique()):
            deme_array = df.loc[df['geography']==j]
            deme_genera_sp = df[['genus', 'species']]
            nb_sp = deme_genera_sp.groupby('genus').nunique().sum().values[0]
            nb_id_max_deme = max(nb_id_max_deme, nb_sp)
    nb_id_max_deme = [nb_id_max_deme] + [0] * nb_time

    # setup stats
    age_ma, msp, vsp, mspd, vspd, mgd, vgd, ng = [], [], [], [], [], [], [], []
    sdt = []
    sdtd = [[] for i in range(nb_deme)]
    nsp = []
    nspd = [[] for i in range(nb_deme)]

    for time in range(nb_time + 1):
        age_ma.append(time)

        df_at_time = df[(df['max_ma'] >= time) & (df['min_ma'] <= time)]  # only the lineages at this time
        if df_at_time.shape[0] != 0:
            # group by genus and species, collects the min and max age of each lineage.
            df_at_time_group = df_at_time.groupby(['genus', 'species']).agg({'max_ma': 'max', 'min_ma': 'min'})
            # Add the median and the difference between min and max.
            df_at_time_group['L'] = df_at_time_group['max_ma'].sub(df_at_time_group['min_ma'], axis=0)
            df_at_time_group['median'] = (df_at_time_group['max_ma'] + df_at_time_group['min_ma']) / 2
            sdt.append(std(df_at_time_group['L'].to_list()))  # Add standard deviation around time
            nsp.append(df_at_time_group['L'].shape[0])  # Add number of species at this time
            ng.append(df_at_time['genus'].nunique())  # add Number of genera at this time

            nb_sp_deme = []
            nb_genera_deme = []

            # do the same as previously but for each deme
            df_at_time['L'] = df_at_time['max_ma'].sub(df_at_time['min_ma'], axis=0)
            for deme in range(nb_deme):
                df_deme_time = df_at_time[df_at_time['geography'] == deme]
                sdtd[deme].append(std(df_deme_time['L'].to_list()))
                nspd[deme].append(df_deme_time.shape[0])
                nb_sp_deme.append(df_deme_time['species'].count())
                nb_genera_deme.append(df_deme_time['genus'].nunique())

            # add the mean and variance for the number of lineages and the number of genera
            mspd.append(mean(nb_sp_deme))
            if len(nb_sp_deme) > 1:
                vspd.append(var(nb_sp_deme, ddof=1))
            else:
                vspd.append(nan)
            mgd.append(mean(nb_genera_deme))
            if len(nb_genera_deme) > 1:
                vgd.append(var(nb_genera_deme, ddof=1))
            else:
                vgd.append(nan)

            nb_sp_genus = []
            df_at_time_simple = df_at_time.drop(
                ['simulation', 'geography', 'specimen', 'extinction_statu', 'max_ma', 'min_ma'], axis=1)
            for genus in df_at_time_simple['genus'].unique():
                df_genus_time = df_at_time_simple[df_at_time_simple['genus'] == genus]
                nb_sp_genus.append(df_genus_time['species'].nunique())
            msp.append(mean(nb_sp_genus))
            if len(nb_sp_genus) > 1:
                vsp.append(var(nb_sp_genus, ddof=1))
            else:
                vsp.append(nan)
        else:  # At this time there are no species
            ng.append(df_at_time['genus'].nunique())
            msp.append(nan)
            mspd.append(nan)
            vsp.append(nan)
            vspd.append(nan)
            mgd.append(nan)
            vgd.append(nan)
            sdt.append(nan)
            nsp.append(nan)
            for deme in range(nb_deme):
                sdtd[deme].append(nan)
                nspd[deme].append(nan)


    # create pandas DataFrame to easily create csv
    df_out = pd.DataFrame({
        'age_Ma' : age_ma,
        'msp' : msp,
        'vsp' : vsp,
        'mspd' : mspd,
        'vspd' : vspd,
        'mgd' : mgd,
        'vgd' : vgd,
        'ng' : ng,
        'sdt' : sdt,
        'sdtd0' : sdtd[0],
        'sdtd1' : sdtd[1],
        'sdtd2' : sdtd[2],
        'sdtd3' : sdtd[3],
        'sdtd4' : sdtd[4],
        'sdtd5' : sdtd[5],
        'sdtd6' : sdtd[6],
        'sdtd7' : sdtd[7],
        'sdtd8' : sdtd[8],
        #'sdtd9' : sdtd[9],
        'nsp' : nsp,
        'nspd0' : nspd[0],
        'nspd1' : nspd[1],
        'nspd2' : nspd[2],
        'nspd3' : nspd[3],
        'nspd4' : nspd[4],
        'nspd5' : nspd[5],
        'nspd6' : nspd[6],
        'nspd7' : nspd[7],
        'nspd8' : nspd[8],
        #'nspd9' : nspd[9],
        'nb_id_max' : nb_id_max,
        'nb_id_max_deme' : nb_id_max_deme[0]
    })

    scenario = os.path.basename(data)
    scenario = scenario.split('_')[2]
    out = fileout+f'/{scenario}_{simul}.csv'
    df_out.to_csv(path_or_buf=out, sep=',', index=False)

def pbdb_to_dataframes(data, nb_time, simul, fileout, nb_deme, sep=','):
    """
    Process the data to create file easily read by R to plot the informations.

    Parameters
    ----------
    data : str
        The path of the file produce by the simulator.
    nb_time : int
        The number of time in simulation.
    simul : int
        The simulation number
    fileout : int
        The path of tje output file produce by this function.
    nb_deme : int
        The number of deme in simulation.
    Returns
    -------
    None
    """
    # read the data
    df = pd.read_csv(data, sep=sep, header=0)

    # Number of unique species during simulation
    df_genera_sp = df[['genus','species']]
    nb_sp = df_genera_sp.groupby('genus').nunique().sum().values[0]
    nb_id_max = [nb_sp] + [0] * nb_time
    nb_id_max_deme = 0

    continent = {}
    count = 0
    # Max number of unique species in a deme during simulation
    for j in df['geography'].unique():
        continent[count] = j
        count += 1
        deme_array = df.loc[df['geography']==j]
        deme_genera_sp = deme_array[['genus', 'species']]
        nb_sp = deme_genera_sp.groupby('genus').nunique().sum().values[0]
        nb_id_max_deme = max(nb_id_max_deme, nb_sp)
    nb_id_max_deme = [nb_id_max_deme] + [0] * nb_time

    print(continent)

    # setup stats
    age_ma, msp, vsp, mspd, vspd, mgd, vgd, ng = [], [], [], [], [], [], [], []
    sdt = []
    sdtd = [[] for i in range(nb_deme)]
    nsp = []
    nspd = [[] for i in range(nb_deme)]

    df[['max_ma', 'min_ma']] = df[['max_ma', 'min_ma']]#.round(0)
    for time in range(nb_time + 1):
        age_ma.append(time)

        df_at_time = df[(df['max_ma'] >= time) & (df['min_ma'] <= time)]  # only the lineages at this time
        if df_at_time.shape[0] != 0:
            # group by genus and species, collects the min and max age of each lineage.
            df_at_time_group = df_at_time.groupby(['genus', 'species']).agg({'max_ma': 'max', 'min_ma': 'min'})
            # Add the median and the difference between min and max.
            df_at_time_group['L'] = df_at_time_group['max_ma'].sub(df_at_time_group['min_ma'], axis=0)
            df_at_time_group['median'] = (df_at_time_group['max_ma'] + df_at_time_group['min_ma']) / 2
            sdt.append(std(df_at_time_group['L'].to_list()))  # Add standard deviation around time
            nsp.append(df_at_time_group.shape[0])  # Add number of species at this time
            ng.append(df_at_time['genus'].nunique())  # add Number of genera at this time

            nb_sp_deme = []
            nb_genera_deme = []

            # do the same as previously but for each deme
            #df_at_time['L'] = df_at_time['max_ma'].sub(df_at_time['min_ma'], axis=0)
            for deme in range(nb_deme):
                if deme in continent:
                    locality = continent[deme]
                    df_deme_time = df_at_time[df_at_time['geography'] == locality]
                    nb_genera_deme.append(df_deme_time['genus'].nunique())
                    df_deme_time = df_deme_time.groupby(['genus', 'species']).agg({'max_ma': 'max', 'min_ma': 'min'})
                    df_deme_time['L'] = df_deme_time['max_ma'].sub(df_deme_time['min_ma'], axis=0)
                    sdtd[deme].append(std(df_deme_time['L'].to_list()))
                    nspd[deme].append(df_deme_time.shape[0])
                    nb_sp_deme.append(df_deme_time.shape[0])
                else:
                    sdtd[deme].append(0)
                    nspd[deme].append(0)
                    nb_sp_deme.append(0)
               

            # add the mean and variance for the number of lineages and the number of genera
            mspd.append(mean(nb_sp_deme))
            if len(nb_sp_deme) > 1:
                vspd.append(var(nb_sp_deme, ddof=1))
            else:
                vspd.append(nan)
            mgd.append(mean(nb_genera_deme))
            if len(nb_genera_deme) > 1:
                vgd.append(var(nb_genera_deme, ddof=1))
            else:
                vgd.append(nan)

            nb_sp_genus = []
            df_at_time_simple = df_at_time.drop(
                ['geography', 'max_ma', 'min_ma'], axis=1)
            for genus in df_at_time_simple['genus'].unique():
                df_genus_time = df_at_time_simple[df_at_time_simple['genus'] == genus]
                nb_sp_genus.append(df_genus_time['species'].nunique())
            msp.append(mean(nb_sp_genus))
            if len(nb_sp_genus) > 1:
                vsp.append(var(nb_sp_genus, ddof=1))
            else:
                vsp.append(nan)
        else:  # At this time there are no species
            ng.append(df_at_time['genus'].nunique())
            msp.append(nan)
            mspd.append(nan)
            vsp.append(nan)
            vspd.append(nan)
            mgd.append(nan)
            vgd.append(nan)
            sdt.append(nan)
            nsp.append(nan)
            for deme in range(nb_deme):
                sdtd[deme].append(nan)
                nspd[deme].append(nan)


    # create pandas DataFrame to easily create csv
    df_out = pd.DataFrame({
        'age_Ma' : age_ma,
        'msp' : msp,
        'vsp' : vsp,
        'mspd' : mspd,
        'vspd' : vspd,
        'mgd' : mgd,
        'vgd' : vgd,
        'ng' : ng,
        'sdt' : sdt,
        'sdtd0' : sdtd[0],
        'sdtd1' : sdtd[1],
        'sdtd2' : sdtd[2],
        'sdtd3' : sdtd[3],
        'sdtd4' : sdtd[4],
        'sdtd5' : sdtd[5],
        'sdtd6' : sdtd[6],
        'sdtd7' : sdtd[7],
        'sdtd8' : sdtd[8],
        'nsp' : nsp,
        'nspd0' : nspd[0],
        'nspd1' : nspd[1],
        'nspd2' : nspd[2],
        'nspd3' : nspd[3],
        'nspd4' : nspd[4],
        'nspd5' : nspd[5],
        'nspd6' : nspd[6],
        'nspd7' : nspd[7],
        'nspd8' : nspd[8],
        'nb_id_max' : nb_id_max,
        'nb_id_max_deme' : nb_id_max_deme
    })

    out = fileout+f'/{simul}.csv'
    df_out.to_csv(path_or_buf=out, sep=',', index=False)

def add_holes(data, nb_time, simul, fileout, nb_deme, prop_holes, global_holes, size_holes):
    # read the data
    df = pd.read_csv(data, sep='\t', header=0)

    # Periods without records
    
    if global_holes:
        nb_holes = choice(nb_time, replace=False, size=int(round(nb_time*prop_holes, 0)))
        for time in nb_holes:
            if size_holes == 10:
                df[(df['max_ma'] <= time+4) & (df['min_ma'] >= time-5)]=None
                df = df.dropna()
                df.loc[df[(df['max_ma'] <= time+4) & ( df['max_ma']>=time-5) & (df['min_ma'] <= time-5)].index, 'max_ma']= time-6
                df.loc[df[(df['max_ma'] >= time+4) & ( df['min_ma']<=time+4) & (df['min_ma'] >= time-5)].index, 'min_ma']= time+5
            elif size_holes == 1:
                df[(df['max_ma'] <= time) & (df['min_ma'] >= time)]=None
                df = df.dropna()
                df.loc[df[(df['max_ma'] <= time) & ( df['max_ma']>=time) & (df['min_ma'] <= time)].index, 'max_ma']= time-1
                df.loc[df[(df['max_ma'] >= time) & ( df['min_ma']<=time) & (df['min_ma'] >= time)].index, 'min_ma']= time+1

    else:
        nb_holes = choice(nb_time, replace=True, size=int(round((nb_deme*nb_time)*prop_holes, 0)))
        select_demes = choice(nb_deme, replace=True, size=int(round((nb_deme*nb_time)*prop_holes, 0)))
        holes_by_demes = {}
        for i in range(len(select_demes)):
            holes_by_demes[nb_holes[i]] = select_demes[i]
        for time in holes_by_demes.keys():
            if size_holes ==10:
                df[(df['max_ma'] <= time+4) & (df['min_ma'] >= time-5) & (df['geography']==holes_by_demes[time])]=None
                df = df.dropna()
                df.loc[df[(df['max_ma'] <= time+4) & ( df['max_ma']>=time-5) & (df['min_ma'] <= time-5) & (df['geography']==holes_by_demes[time])].index, 'max_ma']= time-6
                df.loc[df[(df['max_ma'] >= time+4) & ( df['min_ma']<=time+4) & (df['min_ma'] >= time-5) & (df['geography']==holes_by_demes[time])].index, 'min_ma']= time+5
            elif size_holes ==1:
                df[(df['max_ma'] <= time) & (df['min_ma'] >= time) & (df['geography']==holes_by_demes[time])]=None
                df = df.dropna()
                df.loc[df[(df['max_ma'] <= time) & ( df['max_ma']>=time) & (df['min_ma'] <= time) & (df['geography']==holes_by_demes[time])].index, 'max_ma']= time-1
                df.loc[df[(df['max_ma'] >= time) & ( df['min_ma']<=time) & (df['min_ma'] >= time) & (df['geography']==holes_by_demes[time])].index, 'min_ma']= time+1

    # Number of unique species during simulation
    df_genera_sp = df[['genus','species']]
    nb_sp = df_genera_sp.groupby('genus').nunique().sum().values[0]
    nb_id_max = [nb_sp] + [0] * nb_time
    nb_id_max_deme = 0

    continent = {}
    count = 0
    # Max number of unique species in a deme during simulation
    for j in df['geography'].unique():
        continent[count] = j
        count += 1
        deme_array = df.loc[df['geography']==j]
        deme_genera_sp = deme_array[['genus', 'species']]
        nb_sp = deme_genera_sp.groupby('genus').nunique().sum().values[0]
        nb_id_max_deme = max(nb_id_max_deme, nb_sp)
    nb_id_max_deme = [nb_id_max_deme] + [0] * nb_time

    # setup stats
    age_ma, msp, vsp, mspd, vspd, mgd, vgd, ng = [], [], [], [], [], [], [], []
    sdt = []
    sdtd = [[] for i in range(nb_deme)]
    nsp = []
    nspd = [[] for i in range(nb_deme)]

    df[['max_ma', 'min_ma']] = df[['max_ma', 'min_ma']].round(0)
    for time in range(nb_time + 1):
        age_ma.append(time)

        df_at_time = df[(df['max_ma'] >= time) & (df['min_ma'] <= time)]  # only the lineages at this time
        if df_at_time.shape[0] != 0:
            # group by genus and species, collects the min and max age of each lineage.
            df_at_time_group = df_at_time.groupby(['genus', 'species']).agg({'max_ma': 'max', 'min_ma': 'min'})
            # Add the median and the difference between min and max.
            df_at_time_group['L'] = df_at_time_group['max_ma'].sub(df_at_time_group['min_ma'], axis=0)
            df_at_time_group['median'] = (df_at_time_group['max_ma'] + df_at_time_group['min_ma']) / 2
            sdt.append(std(df_at_time_group['L'].to_list()))  # Add standard deviation around time
            nsp.append(df_at_time_group.shape[0])  # Add number of species at this time
            ng.append(df_at_time['genus'].nunique())  # add Number of genera at this time

            nb_sp_deme = []
            nb_genera_deme = []

            # do the same as previously but for each deme
            #df_at_time['L'] = df_at_time['max_ma'].sub(df_at_time['min_ma'], axis=0)
            for deme in range(nb_deme):
                if deme in continent:
                    locality = continent[deme]
                    df_deme_time = df_at_time[df_at_time['geography'] == locality]
                    nb_genera_deme.append(df_deme_time['genus'].nunique())
                    df_deme_time = df_deme_time.groupby(['genus', 'species']).agg({'max_ma': 'max', 'min_ma': 'min'})
                    df_deme_time['L'] = df_deme_time['max_ma'].sub(df_deme_time['min_ma'], axis=0)
                    sdtd[deme].append(std(df_deme_time['L'].to_list()))
                    nspd[deme].append(df_deme_time.shape[0])
                    nb_sp_deme.append(df_deme_time.shape[0])
                else:
                    sdtd[deme].append(0)
                    nspd[deme].append(0)
                    nb_sp_deme.append(0)
               

            # add the mean and variance for the number of lineages and the number of genera
            mspd.append(mean(nb_sp_deme))
            if len(nb_sp_deme) > 1:
                vspd.append(var(nb_sp_deme, ddof=1))
            else:
                vspd.append(nan)
            mgd.append(mean(nb_genera_deme))
            if len(nb_genera_deme) > 1:
                vgd.append(var(nb_genera_deme, ddof=1))
            else:
                vgd.append(nan)

            nb_sp_genus = []
            df_at_time_simple = df_at_time.drop(
                ['geography', 'max_ma', 'min_ma'], axis=1)
            for genus in df_at_time_simple['genus'].unique():
                df_genus_time = df_at_time_simple[df_at_time_simple['genus'] == genus]
                nb_sp_genus.append(df_genus_time['species'].nunique())
            msp.append(mean(nb_sp_genus))
            if len(nb_sp_genus) > 1:
                vsp.append(var(nb_sp_genus, ddof=1))
            else:
                vsp.append(nan)
        else:  # At this time there are no species
            ng.append(df_at_time['genus'].nunique())
            msp.append(nan)
            mspd.append(nan)
            vsp.append(nan)
            vspd.append(nan)
            mgd.append(nan)
            vgd.append(nan)
            sdt.append(nan)
            nsp.append(nan)
            for deme in range(nb_deme):
                sdtd[deme].append(nan)
                nspd[deme].append(nan)


    # create pandas DataFrame to easily create csv
    df_out = pd.DataFrame({
        'age_Ma' : age_ma,
        'msp' : msp,
        'vsp' : vsp,
        'mspd' : mspd,
        'vspd' : vspd,
        'mgd' : mgd,
        'vgd' : vgd,
        'ng' : ng,
        'sdt' : sdt,
        'sdtd0' : sdtd[0],
        'sdtd1' : sdtd[1],
        'sdtd2' : sdtd[2],
        'sdtd3' : sdtd[3],
        'sdtd4' : sdtd[4],
        'sdtd5' : sdtd[5],
        'sdtd6' : sdtd[6],
        'sdtd7' : sdtd[7],
        'sdtd8' : sdtd[8],
        'nsp' : nsp,
        'nspd0' : nspd[0],
        'nspd1' : nspd[1],
        'nspd2' : nspd[2],
        'nspd3' : nspd[3],
        'nspd4' : nspd[4],
        'nspd5' : nspd[5],
        'nspd6' : nspd[6],
        'nspd7' : nspd[7],
        'nspd8' : nspd[8],
        'nb_id_max' : nb_id_max,
        'nb_id_max_deme' : nb_id_max_deme
    })

    scenario = os.path.basename(data)
    scenario = scenario.split('_')[2]
    if global_holes:
        out = fileout+f'/{scenario}_{simul}_global.csv'
    else:
        out = fileout+f'/{scenario}_{simul}_local.csv'
    df_out.to_csv(path_or_buf=out, sep=',', index=False)

