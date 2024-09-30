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
    'base': 'X',
    'controller': 'o', 
    'miR': 'P',
    'ts3': '^',        
    'ts5': 'v',
    'dual': 's'
}

def get_metadata(path):
    metadata = pd.read_excel(path)

    # Apply colors
    metadata['color'] = metadata['group'].replace(group_palette)
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==2), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==3), 'color'] = colors['red']
    metadata.loc[(metadata['ts_kind']=='NT'), 'color'] = colors['gray']
    metadata.loc[(metadata['ts_num']==4), 'color'] = metadata.loc[(metadata['ts_num']==4), 'color'].apply(get_light_color)
    metadata.loc[(metadata['group']=='miR') & (metadata['miR_loc']=='UTR'), 'color'] = metadata.loc[(metadata['group']=='miR') & (metadata['miR_loc']=='UTR'), 'color'].apply(get_light_color)

    # Apply marker styles
    metadata['markers'] = metadata['group'].replace(group_markers)
    #metadata.loc[(metadata['group']=='controller') & (metadata['ts_kind']=='NT'), 'markers'] = 'X'
    #metadata.loc[(metadata['group']=='marker'), 'markers'] = 'X'

    return metadata

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

def calculate_bins_stats(df, by=['construct','exp','biorep'], stat_list=[sp.stats.gmean, np.std], num_bins=20):

    # Bin by marker quantiles
    df['bin_marker_quantiles'] = df.groupby(by)['marker'].transform(lambda x: pd.qcut(x, q=num_bins, duplicates='drop'))
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