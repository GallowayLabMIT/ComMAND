"""
Simple helper module that encodes both
base style information (like colors) and
helper functions.

This is intended to be imported by Jupyter
notebooks that are used to render each figure.
"""
import matplotlib.pyplot as plt
import rushd as rd
import pandas as pd
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
    'base': 'D',
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
    metadata.loc[(metadata['group']=='controller') & (metadata['ts_kind']=='NT'), 'markers'] = 'X'

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