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

# Sizes are relevant for paper-size figures
font_sizes = {
    'base_size': 8,
    'subpanel_label': 10,
    'poster_base': 16,
    'poster_subpanel': 24,
}

colors = {
    'red': 'crimson', 'orange': 'darkorange', 'yellow': '#ccb804', 'green': 'olive', 'teal': 'teal',
    'blue': '#1650a1', 'purple': 'purple', 'pink': 'hotpink', 'black': 'black', 'gray': 'grey',
}

def get_light_color(color):
    return sns.light_palette(color, 3)[1]

group_palette = {
    'un': colors['black'],          # untransfected / uninfected
    'marker': colors['black'],      # transfection/infection marker only
    'base': colors['black'],        # CDS only, no miR or TS
    'controller': colors['teal'],   # targeting or non-targeting circuit (miR + TS)
    'miR': colors['pink'],          # miR only
    'ts3': colors['purple'],        # TS only, in 3'UTR
    'ts5': colors['purple'],        # TS only, in 5'UTR
}

group_markers = {
    'un': 'H',
    'marker': 'p',
    'base': 'D',
    'controller': 'o', 
    'miR': 'P',
    'ts3': '^',        
    'ts5': 'v',
}

def get_metadata(path):
    metadata = pd.read_excel(path)

    # Apply colors
    metadata['color'] = metadata['group'].replace(group_palette)
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==2), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==3), 'color'] = colors['red']
    metadata.loc[(metadata['ts_kind']=='NT'), 'color'] = colors['gray']
    metadata.loc[(metadata['ts_num']==4), 'color'] = metadata.loc[(metadata['ts_num']==4), 'color'].apply(get_light_color)

    # Apply marker styles
    metadata['markers'] = metadata['group'].replace(group_markers)
    metadata.loc[(metadata['group']=='controller') & (metadata['ts_kind']=='NT'), 'markers'] = 'X'

    return metadata

def rename_multilevel_cols(index):
    if index[1] == '': return index[0]
    else: return index[0] + '_' + index[1]