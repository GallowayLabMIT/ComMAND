"""
Simple helper module that encodes both
base style information (like colors) and
helper functions.

This is intended to be imported by Jupyter
notebooks that are used to render each figure.
"""
import matplotlib.pyplot as plt
import numpy as np
import rushd as rd
import pandas as pd
import scipy as sp
import seaborn as sns

# Sizes are relevant for paper-size figures
font_sizes = {
    'base_size': 8,
    'smaller_size': 7,
    'subpanel_label': 10,
    'poster_base': 16,
    'poster_subpanel': 24,
}

colors = {
    'red': 'crimson', 'orange': 'darkorange', 'yellow': '#ccb804', 'green': '#78AF56', 'teal': 'teal',
    'blue': '#1650a1', 'purple': '#7A1378', 'pink': 'hotpink', 'black': 'black', 'gray': 'grey',
}

def get_light_color(color):
    return sns.light_palette(color, 3)[1]

def get_dark_color(color):
    return sns.dark_palette(color, 3)[1]

group_palette = {
    'un': colors['black'],          # untransfected / uninfected
    'marker': colors['black'],      # transfection/infection marker only
    'base': colors['black'],        # CDS only, no miR or TS
    'controller': colors['teal'],   # targeting or non-targeting circuit (miR + TS)
    'miR': colors['pink'],          # miR only
    'ts3': colors['purple'],        # TS only, in 3'UTR
    'ts5': colors['purple'],        # TS only, in 5'UTR
    'dual': colors['green'],        # part of dual-transcript vector
}

group_markers = {
    'un': 'H',
    'marker': 'p',
    'base': 'D',
    'controller': 'o', 
    'miR': 'P',
    'ts3': '^',        
    'ts5': 'v',
    'dual': 's'
}

def apply_style_tuning(metadata):

    # Create unified color palette for most tuning
    metadata.loc[(metadata['group']=='marker'), 'color'] = 'black'
    metadata.loc[(metadata['group']=='base'), 'color'] = 'black'

    metadata.loc[(metadata['group']=='miR'), 'color'] = colors['purple']
    metadata.loc[(metadata['group']=='miR') & (metadata['miR'].isin(['miR.FF5','miRE.FF5'])), 'color'] = metadata.loc[(metadata['group']=='miR') & (metadata['miR'].isin(['miR.FF5','miRE.FF5'])), 'color'].apply(get_light_color)

    metadata.loc[(metadata['group'].isin(['ts3','ts5'])), 'color'] = colors['blue']
    metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF6x1','FF6x2','FF6x4'])), 'color'] = metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF6x1','FF6x2','FF6x4'])), 'color'].apply(get_dark_color)
    metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF4x1','FF4x2','FF4x4'])), 'color'] = metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF4x1','FF4x2','FF4x4'])), 'color'].apply(get_light_color)
    metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts']=='FF3x1'), 'color'] = metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts']=='FF3x1'), 'color'].apply(get_light_color).apply(get_light_color)

    metadata.loc[(metadata['group']=='controller') & (metadata['miR']=='miR.FF5'), 'color'] = colors['red']
    metadata.loc[(metadata['group']=='controller') & (metadata['miR']=='miRE.FF5'), 'color'] = colors['green']
    metadata.loc[(metadata['group']=='controller') & (metadata['miR']=='miR.FF4'), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['miR']=='miRE.FF4'), 'color'] = colors['teal']

    metadata.loc[(metadata['group']=='controller') & (metadata['ts_kind']=='NT'), 'color'] = colors['gray']
    metadata.loc[(metadata['group']=='controller') & (metadata['ts_kind']=='NT') & (metadata['miR'].isin(['miR.FF4','miRE.FF4'])), 'color'] = metadata.loc[(metadata['group']=='controller') & (metadata['ts_kind']=='NT') & (metadata['miR'].isin(['miR.FF4','miRE.FF4'])), 'color'].apply(get_light_color)

    metadata.loc[((metadata['group']=='controller') & (metadata['ts_num']==2)), 'color'] = metadata.loc[((metadata['group']=='controller') & (metadata['ts_num']==2)), 'color'].apply(get_light_color)
    metadata.loc[((metadata['group']=='controller') & (metadata['ts_num']==4)), 'color'] = metadata.loc[((metadata['group']=='controller') & (metadata['ts_num']==4)), 'color'].apply(get_dark_color)

    metadata.loc[(metadata['ts']=='FF6x3.4x1'), 'color'] = metadata.loc[(metadata['ts']=='FF6x3.4x1'), 'color'].apply(get_light_color)

    # markers
    metadata.loc[(metadata['miR'].isin(['miR.FF5','miR.FF4'])), 'markers'] = 'D'
    metadata.loc[(metadata['miR'].isin(['miRE.FF5','miRE.FF4'])), 'markers'] = 'o'
    metadata.loc[(metadata['group']=='base'), 'markers'] = 'X'

    # linestyles
    metadata['linestyle'] = '-' # default is solid
    metadata.loc[metadata['group']=='marker', 'linestyle'] = ':'
    metadata.loc[(metadata['miR'].isin(['miR.FF5','miR.FF4'])), 'linestyle'] = '--'

    return metadata

def apply_style_promoters(metadata):

    # Create color palette for comparing promoters
    metadata.loc[(metadata['group']=='marker'), 'color'] = 'black'
    metadata.loc[(metadata['group']=='base'), 'color'] = 'black'

    metadata.loc[(metadata['group']=='controller') & (metadata['promoter']=='EF1a'), 'color'] = colors['teal']
    metadata.loc[(metadata['group']=='controller') & (metadata['promoter']=='CAG'), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['promoter']=='EFS'), 'color'] = colors['red']
    metadata.loc[(metadata['group']=='controller') & (metadata['promoter']=='hPGK'), 'color'] = colors['green']

    metadata.loc[(metadata['ts_kind']=='NT'), 'color'] = colors['gray']
    metadata.loc[(metadata['ts_kind']=='NT') & (metadata['promoter']=='EF1a'), 'color'] = metadata.loc[(metadata['ts_kind']=='NT') & (metadata['promoter']=='EF1a'), 'color'].apply(get_light_color)
    metadata.loc[(metadata['ts_kind']=='NT') & (metadata['promoter']=='EFS'), 'color'] = metadata.loc[(metadata['ts_kind']=='NT') & (metadata['promoter']=='EFS'), 'color'].apply(get_dark_color)
    metadata.loc[(metadata['ts_kind']=='NT') & (metadata['promoter']=='hPGK'), 'color'] = metadata.loc[(metadata['ts_kind']=='NT') & (metadata['promoter']=='hPGK'), 'color'].apply(get_dark_color).apply(get_dark_color)

    # markers
    metadata['markers'] = 'o'
    metadata.loc[(metadata['group']=='base'), 'markers'] = 'X'

    # linestyles
    metadata['linestyle'] = '-' # default is solid
    metadata.loc[metadata['group']=='marker', 'linestyle'] = ':'

    return metadata

def get_metadata(path, style='tuning'):

    metadata = pd.read_excel(path)

    metadata['color'] = metadata['group'].replace(group_palette)
    metadata['markers'] = 'o' #metadata['group'].replace(group_markers)
    if style=='tuning': return apply_style_tuning(metadata)
    elif style=='promoters': return apply_style_promoters(metadata)


def gate_data(df, gates):
    df = df.copy()
    exp = df['exp'].values[0] # the same for entire df, assuming df = data.groupby('exp')
    gates_dict = gates.set_index('exp').to_dict('dict') # format: column -> {index: value}
    marker = gates_dict['marker'][exp]
    df['expressing'] = df[marker] > gates_dict[marker][exp]
    df['marker'] = df[marker]
    df['output'] = df[gates_dict['output'][exp]]
    return df

def rename_multilevel_cols(index):
    if index[1] == '': return index[0]
    else: return index[0] + '_' + index[1]

def get_slope(df):
    slope, intercept, r_value, p_value, stderr = sp.stats.linregress(df['bin_marker_quantiles_median_log'], df['output_gmean_log'])
    result = pd.DataFrame(columns=['slope', 'intercept_log', 'r_value', 'p_value', 'stderr'])
    result.loc[len(result.index)] = [slope, intercept, r_value, p_value, stderr]
    return result

def calculate_bins_stats(df, by=['construct','exp','biorep'], stat_list=[sp.stats.gmean, np.std]):

    # Bin by marker quantiles
    df['bin_marker_quantiles'] = df.groupby(by)['marker'].transform(lambda x: pd.qcut(x, q=20, duplicates='drop'))
    quantiles = df.groupby(by+['bin_marker_quantiles'])['marker'].median().rename('bin_marker_quantiles_median').reset_index()
    df_quantiles = df.merge(quantiles, how='left', on=by+['bin_marker_quantiles'])

    # Population stats
    grouped = df_quantiles.groupby(by=by)
    stats = grouped[['marker','output']].agg(stat_list).reset_index().dropna()

    # Rename columns as 'col_stat'
    stats.columns = stats.columns.map(lambda i: rename_multilevel_cols(i))
    stats['count'] = grouped['output'].count().reset_index()['output']

    # Quantile stats & slope
    df_quantiles['bin_marker_quantiles_median_log'] = df_quantiles['bin_marker_quantiles_median'].apply(np.log10)
    df_quantiles['output_log'] = df_quantiles['output'].apply(np.log10)
    grouped = df_quantiles.groupby(by=by+['bin_marker_quantiles_median'])
    stats_quantiles = grouped[['marker','output']].agg(stat_list).reset_index().dropna()

    # Rename columns as 'col_stat'
    stats_quantiles.columns = stats_quantiles.columns.map(lambda i: rename_multilevel_cols(i))
    stats_quantiles['count'] = grouped['output'].count().reset_index()['output']
    stats_quantiles['bin_marker_quantiles_median_log'] = stats_quantiles['bin_marker_quantiles_median'].apply(np.log10)
    stats_quantiles['output_gmean_log'] = stats_quantiles['output_gmean'].apply(np.log10)

    # Compute slope
    fits = stats_quantiles.groupby(by)[stats_quantiles.columns].apply(get_slope).reset_index()
    #display(fits['intercept_log'].max())
    #fits['intercept'] = fits['intercept_log'].apply(lambda x: 10**x)
    
    return df_quantiles, stats, stats_quantiles, fits