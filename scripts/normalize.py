import numpy
import numpy as np
import pandas as pd
from PIL import Image as im

def to_color(val):
    """
    Multiply a value by 255

    Parameters
    ----------
    val : float
        A float between 0 and 1.

    Returns
    -------
    float : Scaled value
    """
    return val*255
def normalization(val, max):
    """
    Normalise a value to be between 0 and 1.

    Parameters
    ----------
    val : float
        Current value to normalise
    max : float
        Maximum value of this stats

    Returns
    -------
    float : Normalised value
    """
    return (val - 0)/(max-0)
def normalize_array(range_file, csv_file, outfile, nb_deme):
    """
    Normalize all the value in input csv file.

    Parameters
    ----------
    range_file : str
        Path to the file with max values of each stats.
    csv_file : str
        Path to input csv file to be standardised.
    outfile : str
        Path to output normalize csv file.
    nb_deme : int
        Number of localities in simulation

    Returns
    -------
    csv file
    """
    ranges = pd.read_csv(range_file, sep='\t', header=0)
    df = pd.read_csv(csv_file, sep='\t', header=0)

    # Recreate the equivalence between the name in csv and the name in the range file.
    equivalence = {'nb_id_max':'nb_id_max', 'nb_id_max_deme':'nb_id_max_deme', 'mean_sp_max_genera':'msp',
                   'var_sp_max_genera':'vsp', 'mean_sp_max_deme':'mspd', 'var_sp_max_deme': 'vspd', 'mean_genus_max_deme': 'mgd',
                   'var_genus_max_deme':'vgd', 'nb_genus_max':'ng', 'standard_deviation_time' :'sdt',
                   'standard_deviation_time_deme' :'sdtd', 'number_species' :'nsp', 'number_species_deme' :'nspd'}
    values = {}
    for i in ranges.keys():
        if i not in ['time_max', 'nb_id_max', 'nb_id_max_deme']:  # Columns which don't need modifications
            if i == 'standard_deviation_time_deme' or i == 'number_species_deme':
                for j in range(nb_deme):
                    name = equivalence[i]+str(j)
                    values[name] = ranges[i][0]
            else:
                values[equivalence[i]] = ranges[i][0]
    for key in values.keys():  # Normalize each columns
        df[key] = df[key].apply(normalization, max=values[key])

    df.to_csv(outfile, sep='\t', header=True, index=False)

def to_image(array_norm, greyout=None, rgbout=None):
    """
    Create pixel images, in greyscale and/or rgba, from normalize csv.

    Parameters
    ----------
    array_norm : str
        Path to normalize csv file.
    greyout : str
        Path to pixels representation in greyscale. By default : None.
        If is None the images not will be produce
    rgbout : str
        Path to pixels representation in rgba. By default : None.
        If is None the images not will be produce

    Returns
    -------
    png file
    """
    df = pd.read_csv(array_norm, sep='\t', header=0)
    greyscale = df.drop(columns=['age_Ma', 'nb_id_max', 'nb_id_max_deme'])
    for key in greyscale.keys():
        greyscale[key] = greyscale[key].apply(to_color)  # change the value into a color
    greyscale = greyscale.round(0)
    greyscale_np = greyscale.to_numpy(dtype=np.uint8, na_value=0)
    if greyout is not None:
        grey = im.fromarray(greyscale_np)
        grey.save(greyout)

    if rgbout is not None:
        R = greyscale.to_numpy(dtype=np.uint8, na_value=0)
        G = greyscale.to_numpy(dtype=np.uint8, na_value=0)
        B = greyscale.to_numpy(dtype=np.uint8, na_value=0)
        A = greyscale.to_numpy(dtype=np.uint8, na_value=0)

        # Each column will have an unique color variation
        for i in range(R.shape[1]):
            if (i+1)%3 == 1:
                G[:,i] = 0
                B[:, i] = 0
            elif (i+1) % 3 == 2:
                R[:, i] = 0
                B[:, i] = 0
            elif (i+1) % 3 == 0:
                R[:, i] = 0
                G[:, i] = 0
            A[:,i] = 255 * ((i+1)/R.shape[1])  # Each stat have different opaque value.
        RGBA = np.stack((R, G, B, A), axis=-1)
        rgb = im.fromarray(RGBA)
        rgb.save(rgbout)

