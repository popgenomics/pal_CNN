import pandas as pd
from numpy import var, mean, nan

def data_to_dataframes(data, file_range, fileout):
    all_range = pd.read_csv(file_range, sep='\t', header=0)
    df = pd.read_csv(data, sep='\t', header=0)
    age_ma, msp, vsp, mspd, vspd, mgd, vgd, ng = [], [], [], [], [], [], [], []
    
    for time in range(all_range['time_max'][0]+1):
        age_ma.append(time)
        df_at_time = df[(df['max_ma'] >= time) & (df['min_ma'] <= time)]
        ng.append(df_at_time['genus'].nunique())
        
        nb_sp_deme = []
        nb_genera_deme = []
        for deme in df_at_time['geography'].unique():
            df_deme_time = df_at_time[df_at_time['geography'] == deme]
            nb_sp_deme.append(df_deme_time['species'].count())
            nb_genera_deme.append(df_deme_time['genus'].nunique())
        mspd.append(mean(nb_sp_deme))
        if len(nb_sp_deme)>1:
            vspd.append(var(nb_sp_deme, ddof=1))
        else:
            vspd.append(nan)
        mgd.append(mean(nb_genera_deme))
        if len(nb_genera_deme)>1:
            vgd.append(var(nb_genera_deme, ddof=1))
        else:
            vgd.append(nan)

        nb_sp_genus = []
        df_at_time_simple = df_at_time.drop(['simulation', 'geography', 'specimen', 'extinction_statu', 'max_ma', 'min_ma'], axis=1)
        for genus in df_at_time_simple['genus'].unique():
            df_genus_time = df_at_time_simple[df_at_time_simple['genus'] == genus]
            nb_sp_genus.append(df_genus_time['species'].nunique())
        msp.append(mean(nb_sp_genus))
        if len(nb_sp_genus)>1:
            vsp.append(var(nb_sp_genus, ddof=1))
        else:
            vsp.append(nan)

        df_out = pd.DataFrame({
            'age_Ma' : age_ma,
            'msp' : msp,
            'vsp' : vsp,
            'mspd' : mspd,
            'vspd' : vspd,
            'mgd' : mgd,
            'vgd' : vgd,
            'ng' : ng
        })

    df_out.to_csv(path_or_buf=fileout, sep='\t', index=False)
