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

# Font sizes for figures
font_sizes = {'base_size': 8, 'smaller_size': 7, 'subpanel_label': 10,}

rc_params = {
    'axes.titlesize': font_sizes['base_size'], 'axes.labelsize': font_sizes['base_size'], 
    'xtick.labelsize': font_sizes['smaller_size'], 'ytick.labelsize': font_sizes['smaller_size'],
    'pdf.fonttype': 42, # exports fonts in pdf
    'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.minor.size': 2, 'ytick.minor.size': 2, 
    'xtick.major.pad': 2, 'ytick.major.pad': 2, 'axes.labelpad': 2,
    'lines.linewidth': 1,
    'axes.spines.right': False, 'axes.spines.top': False, 
}

rc_context = {'font.size': font_sizes['base_size'], 'font.family': 'sans-serif', 'font.sans-serif':['Arial']}

# Figure width guidelines from CellPress, in inches
figure_width = {'full': 6.8504, '1.5-column': 4.48819, '1-column': 3.34646}

# Aesthetics for plotting
#scatter_kwargs = dict(s=4, jitter=0.2, linewidth=0.5, edgecolor='white', legend=False,)
scatter_kwargs = dict(linestyle='none', errorbar='ci', markersize=4, markeredgewidth=0.5, markeredgecolor='white',
                      err_kws=dict(linewidth=1, zorder=0), legend=False,)
line_kwargs = dict(legend=False, dashes=False,  markersize=4, markeredgewidth=0.5, estimator=sp.stats.gmean, 
                   errorbar=lambda x: (sp.stats.gmean(x) / sp.stats.gstd(x), sp.stats.gmean(x) * sp.stats.gstd(x)))
annotate_kwargs = dict(test='t-test_ind', text_format='star', loc='inside', line_height=0, line_width=0.5, text_offset=-2)

# Main color palette
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

def apply_style_designs(metadata):

    # Create color palette for comparing promoters
    metadata.loc[(metadata['group']=='marker'), 'color'] = 'black'
    metadata.loc[(metadata['group']=='base'), 'color'] = 'black'

    metadata.loc[(metadata['group']=='controller') & (metadata['design']==1), 'color'] = colors['teal']
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==2), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==3), 'color'] = colors['red']
    metadata.loc[(metadata['group']=='controller') & (metadata['design']==0), 'color'] = 'black'

    metadata.loc[(metadata['ts_kind']=='NT'), 'color'] = colors['gray']
    metadata.loc[(metadata['ts_kind']=='NT') & (metadata['design']==1), 'color'] = metadata.loc[(metadata['ts_kind']=='NT') & (metadata['design']==1), 'color'].apply(get_light_color)
    metadata.loc[(metadata['ts_kind']=='NT') & (metadata['design']==3), 'color'] = metadata.loc[(metadata['ts_kind']=='NT') & (metadata['design']==3), 'color'].apply(get_dark_color)

    metadata.loc[(metadata['group']=='miR'), 'color'] = colors['purple']
    metadata.loc[(metadata['group']=='miR') & (metadata['miR'].isin(['miR.FF5','miRE.FF5'])), 'color'] = metadata.loc[(metadata['group']=='miR') & (metadata['miR'].isin(['miR.FF5','miRE.FF5'])), 'color'].apply(get_light_color)

    metadata.loc[(metadata['group'].isin(['ts3','ts5'])), 'color'] = colors['blue']
    metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF6x1','FF6x2','FF6x4'])), 'color'] = metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF6x1','FF6x2','FF6x4'])), 'color'].apply(get_dark_color)
    metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF4x1','FF4x2','FF4x4'])), 'color'] = metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts'].isin(['FF4x1','FF4x2','FF4x4'])), 'color'].apply(get_light_color)
    metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts']=='FF3x1'), 'color'] = metadata.loc[(metadata['group'].isin(['ts3','ts5'])) & (metadata['ts']=='FF3x1'), 'color'].apply(get_light_color).apply(get_light_color)

    # markers
    metadata['markers'] = 'o'
    metadata.loc[(metadata['group']=='base'), 'markers'] = 'X'

    # linestyles
    metadata['linestyle'] = '-' # default is solid
    metadata.loc[metadata['group']=='marker', 'linestyle'] = ':'

    return metadata

def apply_style_applications(metadata):

    # Create color palette for comparing promoters & therapeutic genes
    metadata.loc[(metadata['group']=='marker'), 'color'] = 'black'
    metadata.loc[(metadata['group']=='base'), 'color'] = 'black'

    metadata.loc[(metadata['group']=='controller') & (metadata['promoter']=='CAG'), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['name'].str.contains('FMRP')), 'color'] = colors['red']
    metadata.loc[(metadata['group']=='controller') & (metadata['name'].str.contains('FXN')), 'color'] = colors['orange']
    metadata.loc[(metadata['group']=='controller') & (metadata['name'].str.contains('Cre')), 'color'] = colors['green']

    metadata.loc[(metadata['ts_kind']=='NT'), 'color'] = colors['gray']
    metadata.loc[(metadata['ts_kind']=='NT') & (metadata['name'].str.contains('FMRP')), 'color'] = metadata.loc[(metadata['ts_kind']=='NT') & (metadata['name'].str.contains('FMRP')), 'color'].apply(get_light_color)

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
    metadata['ts_label'] = metadata['ts_kind'].replace({'na': 'base', 'NT': 'OL', 'T': 'CL', '–': 'none'})

    if style=='tuning': return apply_style_tuning(metadata)
    elif style=='promoters': return apply_style_promoters(metadata)
    elif style=='designs': return apply_style_designs(metadata)
    elif style=='applications': return apply_style_applications(metadata)


def gate_data(df, gates, group='exp'):
    df = df.copy()
    exp = df[group].values[0] # the same for entire df, assuming df = data.groupby('exp')
    gates_dict = gates.set_index(group).to_dict('dict') # format: column -> {index: value}
    marker = gates_dict['marker'][exp]
    df['expressing'] = df[marker] > gates_dict[marker][exp]
    df['marker'] = df[marker]
    df['output'] = df[gates_dict['output'][exp]]
    return df

def rename_multilevel_cols(index):
    if index[1] == '': return index[0]
    else: return index[0] + '_' + index[1]

def get_slope(df, x='bin_marker_quantiles_median_log', y='output_gmean_log'):
    slope, intercept, r_value, p_value, stderr = sp.stats.linregress(df[x], df[y])
    result = pd.DataFrame(columns=['slope', 'intercept', 'r_value', 'p_value', 'stderr'])
    result.loc[len(result.index)] = [slope, intercept, r_value, p_value, stderr]
    return result

def calculate_bins_stats(df, by=['construct','exp','biorep'], stat_list=[sp.stats.gmean, np.std, sp.stats.variation], num_bins=20):

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
    stats = stats.merge(fits, how='left', on=by)
    
    return df_quantiles, stats

def truncated_poisson(x, mu):
    lower_bound = 1
    # from julia
    #   lcdf = isinf(l) ? 0.0 : cdf(d, l)
    #   ucdf = isinf(u) ? 1.0 : cdf(d, u)
    #   tp = ucdf - lcdf
    #   pdf{T<:Real}(d::Truncated, x::T) = d.lower <= x <= d.upper ? pdf(d.untruncated, x) / d.tp : zero(T)
    truncated_region = 1 - sp.stats.poisson.cdf(lower_bound, mu)
    return sp.stats.poisson.pmf(x, mu) / truncated_region

def map_biorep(df):
    biorep_map = {val:i for i,val in enumerate(df['exp'].unique())}
    d = df.copy()
    d['biorep'] = d['exp'].map(biorep_map)
    return d

# Compute fold changes relative to no-TS conditions
def get_fc(df):
    d = df.copy()
    baseline = d.loc[d['ts']=='none', 'output_gmean'].mean()
    d['output_gmean_fc'] = d['output_gmean'] / baseline
    return d

# Normalize parameter values such that original (middle) value = 1
def normalize_param_val(df):
    d = df.copy()
    vals = d['param_val'].unique()
    d['param_val_norm'] = d['param_val'] / (sorted(vals)[int(len(vals)/2)])
    return d

# Compute instantaneous slope for each param & param_val 
def get_slope_instant(df, x, y):
    d = df.sort_values(x)
    return (list(d[y])[-1] - list(d[y])[-2]) / (list(d[x])[-1] - list(d[x])[-2])

# Get GAPDH Cp
def get_housekeeping_qpcr(df, housekeeping='GAPDH'):
    value = df.loc[df['primers']==housekeeping, 'Cp']
    if value.empty: result = pd.NA
    else: result = value.median()
    df[housekeeping+'_Cp'] = result
    return df

# Get Cp of negative controls (no-cDNA, no-RT)
def get_negative_qpcr(df):
    value = df.loc[df['construct'].isin(['no-RT','no-cDNA']), 'Cp']
    if value.empty: result = 35
    else: result = value.median()
    df['negative_Cp'] = result
    return df

# Get base gene delta_Cp
def get_control_qpcr(df):
    value = df.loc[df['group']=='base', 'delta_Cp']
    if value.empty: result = pd.NA
    else: result = value.values[0]
    df['base_delta_Cp'] = result
    return df

def modify_norm_factor(df):
    param = df['param'].values[0]
    d = df.copy()
    if param in ['α_im', 'α_p']: 
        d['norm_factor'] = d['base_norm_factor'] * d['param_val_norm'].values[0]
    elif param in ['δ_m', 'δ_p']: 
        d['norm_factor'] = d['base_norm_factor'] / d['param_val_norm'].values[0]
    else:
        d['norm_factor'] = d['base_norm_factor']
    return d['norm_factor']

# `load_plates` loads and gates data, renaming channels as needed
def load_plates(data_group, base_path):
    if data_group == 'tuning': return load_plates_tuning(base_path)
    elif data_group == 'ts_num': return load_plates_ts_num(base_path)
    elif data_group == 'miR_characterization': return load_plates_miR_characterization(base_path)
    elif data_group == 'plasmid_titration': return load_plates_plasmid_titration(base_path)
    elif data_group == 'two_gene': return load_plates_two_gene(base_path)
    elif data_group == 'piggybac': return load_plates_piggybac(base_path)
    elif data_group == 'straight-in': return load_plates_straight_in(base_path)
    elif data_group == 'lenti_293T_MEF': return load_plates_lenti_293T_mef(base_path)
    elif data_group == 'lenti_tcell': return load_plates_lenti_tcell(base_path)
    elif data_group == 'lenti_neuron': return load_plates_lenti_neuron(base_path)
    elif data_group == 'lenti_iPS11': return load_plates_lenti_ips11(base_path)
    elif data_group == 'iPS11_transfection': return load_plates_ips11_transfection(base_path)
    elif data_group == 'therapeutic_transfection': return load_plates_therapeutic_transfection(base_path)
    elif data_group == 'lenti_therapeutic': return load_plates_lenti_therapeutic(base_path)
    elif data_group == 'qPCR': return load_plates_qpcr(base_path)
    else: print(f'{data_group} is not a valid data group to load plates.')

def load_plates_tuning(base_path):
    exp90_path = base_path/'kasey'/'2024.03.31_exp90'/'export'
    exp90_2_path = base_path/'kasey'/'2024.04.02_exp90.2'/'export'
    exp90_3_path = base_path/'kasey'/'2024.04.02_exp90.3'/'export'
    exp90_4_path = base_path/'kasey'/'2024.04.05_exp90.4'/'export'
    exp91_path = base_path/'kasey'/'2024.04.08_exp91'/'export'
    exp92_path = base_path/'kasey'/'2024.04.12_exp92'/'export'

    plates = pd.DataFrame({
        'data_path': [exp90_path/'plate1', exp90_path/'plate2', 
                    exp90_2_path, exp90_4_path,
                    exp90_3_path/'plate1', exp90_3_path/'plate2', 
                    exp91_path/'plate1.1', exp91_path/'plate1.2', exp91_path/'plate1.3', 
                    exp91_path/'plate2.1', exp91_path/'plate2.2', exp91_path/'plate2.3',
                    exp92_path/'plate1.1', exp92_path/'plate1.2', exp92_path/'plate1.3', 
                    exp92_path/'plate2.1', exp92_path/'plate2.2', exp92_path/'plate2.3',],
        
        'yaml_path': ([exp90_path/'exp90_plate1_wells.yaml', exp90_path/'exp90_plate2_wells.yaml', 
                    exp90_path/'exp90_plate2_wells.yaml', exp90_path/'exp90_plate1_wells.yaml',
                    exp90_path/'exp90_plate1_wells.yaml', exp90_path/'exp90_plate2_wells.yaml', ] +
                    [exp91_path/'exp91_plate1_wells.yaml']*3 + 
                    [exp91_path/'exp91_plate2.1_wells.yaml', exp91_path/'exp91_plate2.2_wells.yaml', exp91_path/'exp91_plate2.3_wells.yaml'] +
                    [exp92_path/'exp92_plate1_wells.yaml', exp92_path/'exp92_plate1.2_wells.yaml', exp92_path/'exp92_plate1_wells.yaml',
                    exp92_path/'exp92_plate2_wells.yaml', exp92_path/'exp92_plate2.2_wells.yaml', exp92_path/'exp92_plate2_wells.yaml',]
                    ),
        
        'biorep': ([1, 1, 
                    2, 2, 
                    3, 3,] + 
                    [1, 2, 3,]*4),
        
        'exp': (['exp90', 'exp90', 
                'exp90.2', 'exp90.4', 
                'exp90.3', 'exp90.3',] + 
                ['exp91']*6 + 
                ['exp92']*6)
    })
    
    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Add missing gates
    gates.loc[len(gates.index)] = ['exp90.4',0,0,]  
    gates.loc[gates['exp']=='exp90.4', channel_list] = gates.loc[gates['exp']=='exp90.2', channel_list].values

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')
    
    return data

def load_plates_ts_num(base_path):

    exp_path = base_path/'kasey'/'2024.12.04_exp092.3'/'export'

    plates = pd.DataFrame({
        'data_path': [exp_path],
        'yaml_path': [exp_path/'wells.yaml'],
        'biorep': [4],
        'exp': ['exp092.3',]
    })
    
    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')
    
    return data
    
def load_plates_plasmid_titration(base_path): 
    exp100_path = base_path/'kasey'/'2024.07.18_exp100'/'export'

    plates = pd.DataFrame({
        'data_path': [exp100_path, base_path/'kasey'/'2024.09.11_exp100.3'/'export', base_path/'kasey'/'2024.09.28_exp100.5'/'export'],
        'yaml_path': [exp100_path/'wells.yaml']*3,
        'biorep': [1,2,3],
        'exp': ['exp100.1','exp100.3','exp100.5'],
    })

    channel_list = ['mRuby2-A','tagBFP-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Add missing gate for exp100
    gates.loc[len(gates.index)] = ['exp100.1', gates['mGL-A'].mean(), gates['mRuby2-A'].mean(), gates['tagBFP-A'].mean()] 

    # Indicate which channels are relevant for each experiment
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'
    gates['filler'] = 'tagBFP-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')

    return data

def load_plates_miR_characterization(base_path): 
    exp11_path = base_path/'Emma'/'2022.10.11_EXP11'/'Data'
    exp11_controls_path = base_path/'Emma'/'2022.10.11_EXP10'/'data_controls'
    exp49_path = base_path/'Emma'/'2024.04.06_EXP11_replicates'/'Plate_1_EXP49'/'data_singlets'
    exp50_path = base_path/'Emma'/'2024.04.06_EXP11_replicates'/'Plate_2_EXP50'/'data_singlets'
    exp49_50_controls_path = base_path/'Emma'/'2024.04.06_EXP11_replicates'/'Plate_3_Controls'/'data_singlets'

    plates = pd.DataFrame({
        'data_path': [base_path/'Emma'/'2022.10.04_EXP9'/'Data',
                    base_path/'Emma'/'2023.01.16_EXP12'/'Data',
                    base_path/'Emma'/'2023.02.09_EXP13'/'Data',
                    exp11_path, exp11_controls_path, 
                    exp49_path, exp50_path, 
                    exp49_50_controls_path, exp49_50_controls_path],
        
        'yaml_path': ([base_path/'Emma'/'2022.10.04_EXP9'/'Data'/'wells_KL.yaml']*3 + 
                    [exp11_path/'wells_KL.yaml', exp11_controls_path/'wells_KL.yaml', 
                    exp11_path/'wells_KL.yaml', exp11_path/'wells_KL.yaml', 
                    exp49_50_controls_path/'wells_KL.yaml', exp49_50_controls_path/'wells2_KL.yaml', ]),
        
        'biorep': [1,2,3,
                1,1,
                2,3,
                2,3],

        'exp': ['ELP_exp09', 'ELP_exp12', 'ELP_exp13',
                'ELP_exp11', 'ELP_exp11',
                'ELP_exp49', 'ELP_exp50',
                'ELP_exp49', 'ELP_exp50',],
    })
    
    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)
    
    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['ts_construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mRuby2-A'
    gates['output'] = 'mGL-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['ts_construct']!='UT')
    
    return data

def load_plates_two_gene(base_path):
    plates = pd.DataFrame({
        'data_path': [base_path/'kasey'/'2024.04.14_exp93'/'export'/f'plate{n}' for n in range(1,4)] + [base_path/'kasey'/'2024.10.21_exp093.2'/'export'],
        'yaml_path': [base_path/'kasey'/'2024.04.14_exp93'/'export'/'exp93_wells.yaml']*3 + [base_path/'kasey'/'2024.10.21_exp093.2'/'export'/'wells.yaml'],
        'biorep': [1, 2, 3, 4],
        'exp': ['exp093']*3 + ['exp093.2'],
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A','tagBFP-A','SNAP-647-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Rename far-red channel
    data.rename(columns={'SNAP-647-A': 'iRFP670-A'}, inplace=True)
    channel_list = ['mRuby2-A','mGL-A','tagBFP-A','iRFP670-A']
    
    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='GEEC555'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)
    
    # Add manual iRFP670 gate (forgot to include untransfected well)
    gate_iRFP = 2.5e2
    gates['iRFP670-A'] = [gate_iRFP]*len(data['exp'].unique())

    # Indicate which channels are relevant for each experiment (same for both exp093 & exp093.2)
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'iRFP670-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')
    
    return data

def load_plates_piggybac(base_path):
    data_paths = [base_path/'kasey'/'2024.11.07_exp119'/'export', base_path/'kasey'/'2024.11.10_exp119_2'/'export', 
                  base_path/'kasey'/'2024.11.12_exp119.2'/'export', base_path/'kasey'/'2024.11.14_exp119.3'/'export']
    
    plates1 = pd.DataFrame({
        'data_path': data_paths[:2],
        'yaml_path': [p/'wells.yaml' for p in data_paths[:2]],
        'exp': ['exp119', 'exp119_2']
    })
    plates2 = pd.DataFrame({
        'data_path': data_paths[2:],
        'yaml_path': [p/'wells.yaml' for p in data_paths[2:]],
        'exp': ['exp119.2', 'exp119.3']
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data1 = rd.flow.load_groups_with_metadata(plates1, columns=channel_list)
    data2 = rd.flow.load_groups_with_metadata(plates2, columns=channel_list+['mCherry-A'])
    data2.rename(columns={'mRuby2-A': 'YL1-A', 'mCherry-A': 'mRuby2-A'}, inplace=True)
    data = pd.concat([data1,data2], ignore_index=True)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Manually draw gates
    gates = pd.DataFrame({
        'mRuby2-A': [2e2]*4,
        'mGL-A': [2e2]*4,
        'exp': ['exp119', 'exp119_2', 'exp119.2', 'exp119.3']
    })
    
    # Indicate which channels are relevant for each experiment
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing']
    
    return data

def load_plates_straight_in(base_path):

    plates = pd.DataFrame({
        'data_path': [base_path/'Albert'/'Exp14'/'export', base_path/'Albert'/'Exp15'/'export', base_path/'Albert'/'Exp16'/'export',],
        'yaml_path': [base_path/'Albert'/'Exp14'/'export'/'wells.yaml']*3,
        'exp': ['Exp14', 'Exp15', 'Exp16'],
        'biorep': [1,2,3]
    })
    
    # Load data
    channel_list = ['mScarlet2-A','mGreenLantern-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates manually 
    gates = pd.DataFrame({
        'mScarlet2-A': [3e3]*3,
        'mGreenLantern-A': [1e4]*3,
        'exp': ['Exp14', 'Exp15', 'Exp16']
    })

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGreenLantern-A'
    gates['output'] = 'mScarlet2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing']
    
    return data

def load_plates_lenti_293T_mef(base_path):
    base_path_1 = base_path/'kasey'/'2024.04.05_exp89'/'export'
    base_path_2 = base_path/'chris'/'2024.06.02-exp95-lenti-miR-iFFL'/'export'
    plate_list = ['_'.join(x) for x in zip(
            ['plate'+str(i) for i in range(1,10)], 
            (['293T']*3 + ['MEF2A']*3 + ['MEF8A']*3),
            ['P9','P14','P15']*3
    )]

    plates = pd.DataFrame({
        'data_path': [base_path_1/'293T_control', 
                    base_path_1/'293T_plate1', base_path_1/'293T_plate2', base_path_1/'293T_plate3',
                    base_path_1/'MEF_3_plate1', 
                    base_path_1/'MEF_4-1_plate1', base_path_1/'MEF_4-1_plate2', base_path_1/'MEF_4-1_plate3'] +
                    [base_path_2/p for p in plate_list],
        'yaml_path': [base_path_1/'kasey_yaml2'/'plate_control.yaml', 
                    base_path_1/'kasey_yaml2'/'plate01.yaml', base_path_1/'kasey_yaml2'/'plate02.yaml', base_path_1/'kasey_yaml2'/'plate03.yaml',
                    base_path_1/'kasey_yaml2'/'mef_3_plate01.yaml', 
                    base_path_1/'kasey_yaml2'/'mef_4-1_plate01.yaml', base_path_1/'kasey_yaml2'/'mef_4-1_plate02.yaml', base_path_1/'kasey_yaml2'/'mef_4-1_plate03.yaml'] +
                    [base_path_2/(p+'_metadata.yaml') for p in plate_list],
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Add more metadata
    data['exp'] = data['cell_type'] + '_' + data['virus_batch']
    data['cell'] = data['cell_type'].apply(lambda x: x.split('-')[0])
    data = data.groupby('cell')[data.columns].apply(map_biorep).reset_index(drop=True)
    data['moi'] = data['virus_dilution']

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[(data['virus_dilution']==0)].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Define gates based on control plate
    gates.loc[len(gates.index)] = ['293T_P10'] + list(gates.loc[gates['exp']=='293T_na', channel_list].mean().values)
    gates.loc[len(gates.index)] = ['293T_P14_'] + list(gates.loc[gates['exp']=='293T_na', channel_list].mean().values)
    gates.loc[len(gates.index)] = ['293T_P16'] + list(gates.loc[gates['exp']=='293T_na', channel_list].mean().values) 
    gates.loc[len(gates.index)] = ['MEF-3_P10'] + list(gates.loc[gates['exp'].str.contains('MEF'), channel_list].mean().values)
    gates.loc[len(gates.index)] = ['MEF-4-1_P10'] + list(gates.loc[gates['exp'].str.contains('MEF'), channel_list].mean().values)
    gates.loc[len(gates.index)] = ['MEF-4-1_P14'] + list(gates.loc[gates['exp'].str.contains('MEF'), channel_list].mean().values)
    gates.loc[len(gates.index)] = ['MEF-4-1_P16'] + list(gates.loc[gates['exp'].str.contains('MEF'), channel_list].mean().values)

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['virus_dilution']!=0)
    
    return data

# ignore base_path for just this one right now, will change later
def load_plates_lenti_tcell(base_path):
    tcell_path = rd.datadir/'instruments'/'data'/'collaborators'/'birnbaum_steph'

    plates = pd.DataFrame({
        'data_path': [tcell_path/'2024-06-10 Galloway Exp 1'/'export', tcell_path/'2024-10-25 Galloway 2'/'export',
                      tcell_path/'2024-11-06 Galloway 3'/'export'],
        'yaml_path': [tcell_path/'2024-06-10 Galloway Exp 1'/'metadata.yaml', tcell_path/'2024-10-25 Galloway 2'/'export'/'metadata.yaml',
                      tcell_path/'2024-11-06 Galloway 3'/'export'/'wells.yaml'],
        'exp': ['steph1', 'steph2', 'steph3']
    })
    
    # Load data
    channel_list = ['FITC-A', 'PE-A', 'APC-A750-A', 'PB450-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Rename channels
    d1 = data[data['exp']=='steph1'].copy()
    d2 = data[data['exp']=='steph2'].copy()
    d3 = data[data['exp']=='steph3'].copy()
    d1 = d1.rename({'FITC-A': 'mGL-A', 'PE-A': 'mRuby2-A', 'APC-A750-A': 'livedead-A'}, axis=1)
    d2 = d2.rename({'FITC-A': 'mGL-A', 'PE-A': 'mRuby2-A', 'PB450-A': 'livedead-A'}, axis=1)
    d3 = d3.rename({'FITC-A': 'mGL-A', 'PE-A': 'mRuby2-A', 'APC-A750-A': 'livedead-A'}, axis=1)
    d1['biorep'] = 1
    d2['biorep'] = 2

    data = pd.concat([d1, d2, d3], ignore_index=True)

    # Remove negative channel values
    channel_list = ['mGL-A', 'mRuby2-A', 'livedead-A']
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp','biorep'])[channel].apply(lambda x: x.quantile(0.9999))
    gates.reset_index(inplace=True)

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['biorep'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Manually adjust marker gates to better isolate expressing population
    gates['mGL-A'] = [4e4, 4e4, 2e5, 2e5]

    # Gate data by marker expression
    data = data.groupby('biorep')[data.columns].apply(lambda x: gate_data(x,gates,'biorep'))
    data.reset_index(inplace=True, drop=True)

    # Gate live cells (livedead-A < manual gate)
    data.loc[data['biorep']==1, 'live'] = data.loc[data['biorep']==1, 'livedead-A'] < 3e3
    data.loc[data['biorep']==2, 'live'] = data.loc[data['biorep']==2, 'livedead-A'] < 5e4
    data.loc[data['biorep']>2, 'live'] = data.loc[data['biorep']>2, 'livedead-A'] < 1e4

    data['gated'] = data['expressing'] & data['live'] & (data['construct']!='UT')
    data['cell'] = 'tcell'
    
    return data

def load_plates_lenti_neuron(base_path):
    neuron_path = base_path/'chris'/'2024.06.15-rat-neurons'

    plates = pd.DataFrame({
        'data_path': [base_path/'chris'/'2024.06.15-rat-neurons'/'export', base_path/'kasey'/'2024.11.12_exp098.2'/'export', base_path/'kasey'/'2024.11.23_exp098.4'/'export'],
        'yaml_path': [base_path/'chris'/'2024.06.15-rat-neurons'/'metadata.yaml', base_path/'kasey'/'2024.11.12_exp098.2'/'export'/'wells.yaml', base_path/'kasey'/'2024.11.23_exp098.4'/'export'/'wells.yaml'],
        'exp': ['exp098', 'exp098.2', 'exp098.4'],
        'cell': ['neuron']*3,
        'dox': [1000]*3
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Adjust marker gate to better isolate infected population
    gates['mGL-A'] = 1e3

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')
    
    return data

def load_plates_lenti_ips11(base_path):
    plates = pd.DataFrame({
        'data_path': [base_path/'kasey'/'2024.10.03_exp116'/'export'/'plate2', base_path/'kasey'/'2024.10.21_exp116.2'/'export'],
        'yaml_path': [base_path/'kasey'/'2024.10.03_exp116'/'export'/'plate2'/'exp116_wells.yaml', base_path/'kasey'/'2024.10.21_exp116.2'/'export'/'wells.yaml'],
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)
    data = data[~data['construct'].isna()]

    # Add more metadata
    data['exp'] = 'exp116_' + data['biorep'].astype(str)
    data['cell'] = 'iPS11'
    data['moi'] = 1

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UI'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Adjust marker gate to better isolate infected population
    gates['mGL-A'] = [2e3]*3

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UI')
    
    return data

def load_plates_ips11_transfection(base_path):
    exp83_5_path = base_path/'kasey'/'2024.04.29_exp83.5'/'export'
    plates = pd.DataFrame({
        'data_path': [exp83_5_path, 
                      base_path/'Emma'/'2024.06.14_EXP65'/'data_singlets', base_path/'Emma'/'2024.06.21_EXP65.2'/'data_singlets'],
        'yaml_path': [exp83_5_path/'exp83.5_wells.yaml']*3,
        'biorep': [1, 2, 3],
        'exp': ['exp83.5', 'elp_exp65', 'elp_exp65.2']
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    d = rd.flow.load_groups_with_metadata(plates, columns=channel_list+['YL1-A'])
    
    # Rename channels
    d1 = d[d['exp']=='exp83.5'].copy()
    d2 = d[d['exp']!='exp83.5'].copy()
    d2.rename(columns={'YL1-A': 'mRuby2-A', 'mRuby2-A': 'mCherry-A'}, inplace=True)
    data = pd.concat([d1, d2], ignore_index=True)
    
    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]
    
    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'
    
    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')
    
    return data

def load_plates_therapeutic_transfection(base_path):
    rep1_2_path = base_path/'Emma'/'2024.06.05_EXP56'/'data_singlets'
    rep3_path = base_path/'Emma'/'2024.06.09_EXP60'/'data_singlets'

    plates = pd.DataFrame({
        'data_path': [rep1_2_path, rep1_2_path, rep3_path],
        'yaml_path': [rep1_2_path/'elp_exp56_biorep_1_wells.yaml', rep1_2_path/'elp_exp56_biorep_2_wells.yaml', rep3_path/'elp_exp56_biorep_3_wells.yaml'],
        'exp': ['elp_exp56.1','elp_exp56.2','elp_exp60'],
        'biorep': [1,2,3]
    })

    # Load data
    channel_list = ['mRuby2-A','EGFP-A','iRFP-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UT'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'iRFP-A'
    gates['output'] = 'mRuby2-A' # only for FXN

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UT')
    
    return data

def load_plates_lenti_therapeutic(base_path):
    exp_path = base_path/'kasey'/'2024.11.18_exp117.5'/'export'
    plates = pd.DataFrame({
        'data_path': [exp_path/'plate1', exp_path/'plate2',],
        'yaml_path': [exp_path/'plate1'/'wells.yaml', exp_path/'plate2'/'wells.yaml'],
        'exp': ['exp117.5']*2
    })

    # Load data
    channel_list = ['mRuby2-A','mGL-A']
    data = rd.flow.load_groups_with_metadata(plates, columns=channel_list)

    # Add more metadata
    data['cell'] = 'therapeutic'
    data['moi'] = 1
    data['marker'] = data['mGL-A']
    data['output'] = data['mRuby2-A']

    # Change output/marker channel labels for FMRP
    fmrp_list = ['RC262','RC263','RC264', 'RC284','RC285','RC286']
    data.loc[data['construct'].isin(fmrp_list), 'marker'] = data.loc[data['construct'].isin(fmrp_list), 'mRuby2-A']
    data.loc[data['construct'].isin(fmrp_list), 'output'] = data.loc[data['construct'].isin(fmrp_list), 'mGL-A']
    data.loc[data['construct'].isin(fmrp_list), 'exp'] = data.loc[data['construct'].isin(fmrp_list), 'exp'] + '_FMRP'

    # Remove negative channel values
    for c in channel_list: data = data[data[c]>0]

    # Draw gates
    gates = pd.DataFrame()
    for channel in channel_list:
        gates[channel] = data[data['construct']=='UI'].groupby(['exp'])[channel].apply(lambda x: x.quantile(0.999))
    gates.reset_index(inplace=True)
    
    # Duplicate gates for artificial FMRP "experiment"
    gates.loc[len(gates.index)] = ['exp117.5_FMRP', gates['mGL-A'].mean(), gates['mRuby2-A'].mean(),] 

    # Adjust marker gate to better isolate infected population
    gates['mGL-A'] = [2e3]*2

    # Indicate which channels are relevant for each experiment
    gates.sort_values(['exp'], inplace=True)
    gates['marker'] = 'mGL-A'
    gates['output'] = 'mRuby2-A'
    gates.loc[gates['exp'].str.contains('FMRP'), 'marker'] = 'mRuby2-A'
    gates.loc[gates['exp'].str.contains('FMRP'), 'output'] = 'mGL-A'

    # Gate data by marker expression
    data = data.groupby('exp')[data.columns].apply(lambda x: gate_data(x,gates))
    data.reset_index(inplace=True, drop=True)
    data['gated'] = data['expressing'] & (data['construct']!='UI')
    
    return data

# ignore base_path for just this one right now, will change later
def load_plates_qpcr(base_path):
    base_path = rd.datadir/'instruments'/'data'/'qPCR'/'emma'/'command'
    plates = pd.DataFrame({
        'data_path': [base_path/'2024.11.13_command'/'2024.11.13_ELP_command_qPCR_Cp.txt', base_path/'2024.11.22_command'/'2024.11.22_ELP_command_qPCR_Cp.txt', base_path/'2024.12.03_command'/'2024.12.03_ELP_command_qPCR_Cp.txt'],
        'yaml_path': [base_path/'2024.11.13_command'/'wells.yaml', base_path/'2024.11.22_command'/'wells.yaml', base_path/'2024.12.03_command'/'wells.yaml'],
        'exp': ['exp123', 'exp123.2', 'exp123.4']
    })

    group_list = []
    for group in plates.to_dict(orient="index").values():

        # Load data in group
        group_data = pd.read_csv(group['data_path'], sep='\t', header=1, usecols=['Pos','Cp'])
        group_data.rename(columns={'Pos': 'well'}, inplace=True)

        # Add associated metadata from `plates` (not paths)
        for k, v in group.items():
            if not (k == "data_path") and not (k == "yaml_path"):
                group_data[k] = v

        # Add well metadata
        metadata = pd.DataFrame.from_dict(rd.flow.load_well_metadata(group['yaml_path']))
        metadata.reset_index(names='well', inplace=True)
        group_data = group_data.merge(metadata, how='left', on='well',)

        group_list.append(group_data)

    # Concatenate all the data into a single DataFrame
    data = pd.concat(group_list, ignore_index=True).replace(np.NaN, pd.NA)
    data.dropna(subset='Cp', inplace=True)

    return data

def load_data(base_path, metadata_path, which, metadata_style='tuning'):
    if which == 'tuning': return load_data_tuning(base_path, metadata_path, metadata_style)
    elif which == 'plasmid_titration': return load_data_plasmid_titration(base_path, metadata_path)
    elif which == 'miR_characterization': return load_data_miR_characterization(base_path, metadata_path)
    elif which == 'two_gene': return load_data_two_gene(base_path, metadata_path)
    elif which == 'piggybac': return load_data_piggybac(base_path, metadata_path) # PiggyBac and STRAIGHT-IN integrations
    elif which == 'lenti': return load_data_lenti(base_path, metadata_path) # all lentivirus data (all cell types + therapeutic genes)
    elif which == 'application': return load_data_application(base_path, metadata_path) # iPS11 and therapeutic gene transfections
    elif which == 'qPCR': return load_data_qpcr(base_path, metadata_path)
    else: print(f'{which} is not a valid data group.')

def load_data_tuning(base_path, metadata_path, metadata_style):

    data_groups = ['tuning', 'ts_num']

    data_list = []
    for data_group in data_groups:
        cache_path = rd.rootdir/'data'/(data_group+'.gzip')
        data = pd.DataFrame()
        if cache_path.is_file(): data = pd.read_parquet(cache_path)
        else: 
            data = load_plates(data_group, base_path)
            data.to_parquet(rd.outfile(cache_path))
        data_list.append(data)

    # Combine data into a single dataframe
    data = pd.concat(data_list, ignore_index=True)

    # Bin data and calculate statistics
    df_quantiles, df_stats = calculate_bins_stats(data[data['gated']])

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx', style=metadata_style)
    data = data.merge(metadata, how='left', on='construct')
    df_quantiles = df_quantiles.merge(metadata, how='left', on='construct')
    df_stats = df_stats.merge(metadata, how='left', on='construct')

    return data, df_quantiles, df_stats, metadata

def load_data_plasmid_titration(base_path, metadata_path):

    # Load and gate raw data
    cache_path = rd.rootdir/'data'/('plasmid_titration.gzip')
    data = pd.DataFrame()
    if cache_path.is_file(): data = pd.read_parquet(cache_path)
    else: 
        data = load_plates('plasmid_titration', base_path)
        data.to_parquet(rd.outfile(cache_path))

    # Bin data and calculate statistics
    df_quantiles, df_stats = calculate_bins_stats(data[data['gated']], by=['construct','amount','exp','biorep'], num_bins=10)

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx')
    data = data.merge(metadata, how='left', on='construct')
    df_quantiles = df_quantiles.merge(metadata, how='left', on='construct')
    df_stats = df_stats.merge(metadata, how='left', on='construct')

    return data, df_quantiles, df_stats, metadata

def load_data_miR_characterization(base_path, metadata_path):

    # Load and gate raw data
    cache_path = rd.rootdir/'data'/('miR_characterization.gzip')
    data = pd.DataFrame()
    if cache_path.is_file(): data = pd.read_parquet(cache_path)
    else: 
        data = load_plates('miR_characterization', base_path)
        data.to_parquet(rd.outfile(cache_path))
    
    data['condition'] = data['miR_construct'] + '_' + data['ts_construct']

    # Bin data and calculate statistics
    df_quantiles, df_stats = calculate_bins_stats(data[data['gated']], by=['condition','miR_construct','ts_construct','biorep','exp'])

    # Add metadata for miR_construct
    metadata1 = pd.read_excel(metadata_path/'miR-metadata.xlsx')
    data = data.merge(metadata1, how='left', on='miR_construct')
    df_quantiles = df_quantiles.merge(metadata1, how='left', on='miR_construct')
    df_stats = df_stats.merge(metadata1, how='left', on='miR_construct')

    # Add metadata for ts_construct
    metadata2 = pd.read_excel(metadata_path/'ts-metadata.xlsx')
    data = data.merge(metadata2, how='left', on='ts_construct')
    df_quantiles = df_quantiles.merge(metadata2, how='left', on='ts_construct')
    df_stats = df_stats.merge(metadata2, how='left', on='ts_construct')

    # Combine metadata
    metadata = data.drop_duplicates('condition')[['miR_construct','ts_construct','condition']]
    metadata = metadata.merge(metadata1, how='left', on='miR_construct')
    metadata = metadata.merge(metadata2, how='left', on='ts_construct')

    # For this experiment, compute fold-change expression of each ts relative to ts=none for each miR
    df_stats = df_stats.groupby(by=['miR_construct','biorep','exp'])[df_stats.columns].apply(get_fc).reset_index(drop=True)

    return data, df_quantiles, df_stats, metadata

def load_data_two_gene(base_path, metadata_path):

    # Load and gate raw data
    cache_path = rd.rootdir/'data'/('two_gene.gzip')
    data = pd.DataFrame()
    if cache_path.is_file(): data = pd.read_parquet(cache_path)
    else: 
        data = load_plates('two_gene', base_path)
        data.to_parquet(rd.outfile(cache_path))

    data['condition'] = data['construct'] + '_' + data['construct2']

    # Bin data and calculate statistics
    df_quantiles, df_stats = calculate_bins_stats(data[data['gated']], by=['condition','construct','construct2','biorep','exp'])

    # Add metadata for construct
    metadata1 = get_metadata(metadata_path/'construct-metadata.xlsx', 'designs')
    data = data.merge(metadata1, how='left', on='construct')
    df_quantiles = df_quantiles.merge(metadata1, how='left', on='construct')
    df_stats = df_stats.merge(metadata1, how='left', on='construct')

    # Add metadata for construct2
    metadata2 = pd.read_excel(metadata_path/'construct2-metadata.xlsx')
    data = data.merge(metadata2, how='left', on='construct2')
    df_quantiles = df_quantiles.merge(metadata2, how='left', on='construct2')
    df_stats = df_stats.merge(metadata2, how='left', on='construct2')

    # Combine metadata
    metadata = data.drop_duplicates('condition')[['construct','construct2','condition']]
    metadata.dropna(inplace=True)
    metadata = metadata.merge(metadata1, how='left', on='construct')
    metadata = metadata.merge(metadata2, how='left', on='construct2')

    # Create color/marker palettes
    metadata.loc[(metadata['gene']=='1T') & (metadata['design']==1), 'color'] = colors['teal']
    metadata.loc[metadata['gene']=='2T', 'color'] = colors['green']
    metadata.loc[metadata['gene']=='2V', 'color'] = colors['purple']

    # markers
    metadata['markers'] = 'X'
    metadata.loc[metadata['gene']=='1T', 'markers'] = 'o'
    metadata.loc[metadata['gene']=='2T', 'markers'] = 'D'
    metadata.loc[metadata['gene']=='2V', 'markers'] = 's'

    ts_label = {'na': 'base', 'NT': 'OL', 'T': 'CL', 'none': '–'}
    metadata['ts_label'] = metadata['ts_kind'].replace(ts_label)
    metadata['kind'] = metadata['gene'] + '_' + metadata['design'].astype(int).astype(str)

    return data, df_quantiles, df_stats, metadata

def load_data_piggybac(base_path, metadata_path):

    data_groups = ['piggybac', 'straight-in']

    data_list = []
    for data_group in data_groups:
        cache_path = rd.rootdir/'data'/(data_group+'.gzip')
        data = pd.DataFrame()
        if cache_path.is_file(): data = pd.read_parquet(cache_path)
        else: 
            data = load_plates(data_group, base_path)
            data.to_parquet(rd.outfile(cache_path))
        data['vector'] = data_group
        data_list.append(data)

    # Combine data into a single dataframe
    data = pd.concat(data_list, ignore_index=True)

    # Bin data and calculate statistics
    df_quantiles, df_stats = calculate_bins_stats(data[data['gated']], by=['construct','exp','biorep','vector'])

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx')
    data = data.merge(metadata, how='left', on='construct')
    df_quantiles = df_quantiles.merge(metadata, how='left', on='construct')
    df_stats = df_stats.merge(metadata, how='left', on='construct')

    return data, df_quantiles, df_stats, metadata

# Combine data from all lenti experiments
def load_data_lenti(base_path, metadata_path):

    data_groups = ['lenti_293T_MEF', 'lenti_tcell', 'lenti_neuron', 'lenti_iPS11', 'lenti_therapeutic']

    data_list = []
    for data_group in data_groups:
        cache_path = rd.rootdir/'data'/(data_group+'.gzip')
        data = pd.DataFrame()
        if cache_path.is_file(): data = pd.read_parquet(cache_path)
        else: 
            data = load_plates(data_group, base_path)
            data.to_parquet(rd.outfile(cache_path))
        data_list.append(data)
    
    # Combine data into a single dataframe
    data = pd.concat(data_list, ignore_index=True)
    
    # Bin data and calculate statistics
    fewer_bins = ['neuron','iPS11']
    df_quantiles, df_stats = calculate_bins_stats(data[(data['gated']) & ~(data['cell'].isin(fewer_bins))].copy(), by=['construct','moi','dox','cell','biorep','exp'])
    df_quantiles2, df_stats2 = calculate_bins_stats(data[(data['gated']) & (data['cell'].isin(fewer_bins))].copy(), by=['construct','moi','dox','cell','biorep','exp'], num_bins=10)

    quantiles = pd.concat([df_quantiles, df_quantiles2], ignore_index=True)
    stats = pd.concat([df_stats, df_stats2], ignore_index=True)

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx', 'designs')
    data = data.merge(metadata, how='left', on='construct')
    quantiles = quantiles.merge(metadata, how='left', on='construct')
    stats = stats.merge(metadata, how='left', on='construct')

    return data, quantiles, stats, metadata

# Combine data from application-relevant transfections
def load_data_application(base_path, metadata_path):

    data_groups = ['therapeutic_transfection', 'iPS11_transfection']

    data_list = []
    for data_group in data_groups:
        cache_path = rd.rootdir/'data'/(data_group+'.gzip')
        data = pd.DataFrame()
        if cache_path.is_file(): data = pd.read_parquet(cache_path)
        else: 
            data = load_plates(data_group, base_path)
            data['data_group'] = data_group
            data.to_parquet(rd.outfile(cache_path))
        data_list.append(data)
    
    # Combine data into a single dataframe
    data = pd.concat(data_list, ignore_index=True)

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx', 'applications')
    data = data.merge(metadata, how='left', on='construct')
    data = data[~data['construct'].isna()]
    
    # Rename output channel for FMRP
    data.loc[data['name'].str.contains('FMRP'), 'output'] = data.loc[data['name'].str.contains('FMRP'), 'EGFP-A']

    # Bin data and calculate statistics
    df_quantiles, df_stats = calculate_bins_stats(data[data['gated']], by=['data_group','construct','biorep','exp'])
    df_stats = df_stats.merge(metadata, how='left', on='construct')

    return data, df_quantiles, df_stats, metadata

def load_data_qpcr(base_path, metadata_path):

    # Load raw (Cp) data
    cache_path = rd.rootdir/'data'/'qPCR.gzip'
    data = pd.DataFrame()
    if cache_path.is_file(): data = pd.read_parquet(cache_path)
    else: 
        data = load_plates('qPCR', base_path)
        data.to_parquet(rd.outfile(cache_path))

    # Remove points (biorep + construct + primers + technical rep) where Cp = 35 
    # (poor amplification, due technical outliers or zero expression conditions)
    # Also, ignore -dox conditions
    print('Remove technical reps with Cp=35')
    display(data[(data['Cp']>=35)])
    data = data[(data['Cp']<35) & (data['dox'])].copy()

    # Calculate GAPDH Cp for each condition (biorep + construct)
    data = data.groupby(['exp','biorep','construct'])[data.columns].apply(get_housekeeping_qpcr).reset_index(drop=True)

    # Combine technical reps (median)
    stats = data.groupby(['exp','biorep','construct','primers','GAPDH_Cp'])[['Cp']].median().reset_index().dropna(subset='Cp')

    # Exclude points (biorep/exp + primers) where Cp ~ Cp of negative controls (no-cDNA, no-RT)
    stats = stats.groupby(['exp','primers'])[stats.columns].apply(get_negative_qpcr).reset_index(drop=True)
    print('Remove points where Cp ~ Cp of negative controls')
    display(stats[(stats['Cp'] >= stats['negative_Cp'])])
    filtered = stats[(stats['Cp'] < stats['negative_Cp']) & ~(stats['construct'].isin(['no-RT','no-cDNA']))]

    # Exclude conditions (biorep + construct) with low GAPDH expression (Cp >= 21)
    print('Remove conditions with low GAPDH expression')
    display(filtered[filtered['GAPDH_Cp'] >= 21])
    filtered = filtered[filtered['GAPDH_Cp'] < 21]
    filtered['expression'] = 2**(-filtered['Cp'])

    # Normalize expression relative to GAPDH for each condition (biorep + construct)
    filtered['delta_Cp'] = filtered['Cp'] - filtered['GAPDH_Cp']
    filtered['norm_expression'] = 2**(-filtered['delta_Cp'])

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx', 'applications')
    data = data.merge(metadata, how='left', on='construct')
    filtered = filtered.merge(metadata, how='left', on='construct')

    # Calculate base delta_Cp for each condition (biorep + primers)
    filtered = filtered.groupby(['exp','biorep','primers'])[filtered.columns].apply(get_control_qpcr).reset_index(drop=True)

    # Exclude conditions (biorep + construct) missing their corresponding base gene 
    # e.g., base gene previously excluded due to low amplification
    print('Remove conditions missing their corresponding base gene')
    display(filtered[filtered['base_delta_Cp'].isna()])
    filtered.dropna(subset='base_delta_Cp', inplace=True)

    # Normalize expression relative to base gene (delta_delta_Cp) for each condition (biorep + primers)
    filtered['delta_delta_Cp'] = filtered['delta_Cp'] - filtered['base_delta_Cp']
    filtered['relative_expression'] = 2**(-filtered['delta_delta_Cp'])
    
    return data, filtered, metadata

# Load modeling results
def load_modeling(base_path, which):
    if which == 'param_sweeps': return load_modeling_param_sweeps(base_path)
    elif which == 'stochastic_sims': return load_modeling_stochastic_sims(base_path)
    elif which == 'missplicing': return load_modeling_missplicing(base_path)
    else: print(f'{which} is not a valid modeling type.')

# later: move modeling results to different location
#   current base_path = rd.rootdir/'output'
def load_modeling_param_sweeps(base_path):

    # Load parameter sweeps from ODE model
    simulation_path = base_path/'modeling'/'julia_param_sweeps'/'per_param'/'sweep_df.gzip'
    data = pd.DataFrame()
    if simulation_path.is_file(): 
        data = pd.read_parquet(simulation_path)

    # Normalize parameter values such that original (middle) value = 1
    data = data.groupby('param')[data.columns].apply(normalize_param_val).reset_index(drop=True)

    # Compute value of unregulated gene
    alpha_rna = 4.67e-2     # params from `miR_iFFL.jl`
    delta_mrna = 2.88e-4
    alpha_p = 3.33e-4
    delta_p = 9.67e-5
    data['unreg'] = data['copy_num'] * (alpha_rna * alpha_p) / (delta_mrna * delta_p)

    # Compute instantaneous slope for each param & param_val 
    col_list = ['copy_num','protein']
    slopes = data.groupby(['param','param_val','param_val_norm'])[data.columns].apply(lambda x: get_slope_instant(x, *col_list)).rename('slope').reset_index()
    slopes['base_norm_factor'] = (delta_mrna * delta_p) / (alpha_rna * alpha_p)

    result = slopes.groupby(['param','param_val_norm'])[slopes.columns].apply(modify_norm_factor).rename('norm_factor').reset_index()
    slopes['norm_factor'] = result['norm_factor']
    slopes['slope_norm'] = slopes['slope'] * slopes['norm_factor']
    slopes['gene'] = 'output'

    # Compute slope for unregulated gene
    by = ['param','param_val','param_val_norm']
    slopes_unreg = data.groupby(by)[data.columns].apply(lambda x: get_slope_instant(x, *['copy_num','unreg'])).rename('slope').reset_index()
    slopes_unreg['norm_factor'] = result['norm_factor']
    slopes_unreg['slope_norm'] = slopes_unreg['slope'] * slopes_unreg['norm_factor']
    slopes_unreg['gene'] = 'marker'

    # Combine into a single dataframe
    slopes = pd.concat([slopes, slopes_unreg])

    return data, slopes

# later: move modeling results to a different location
#    current base_path = rd.datadir/'projects'/'miR-iFFL'
def load_modeling_stochastic_sims(base_path):

    # Load stochastic simulations
    simulation_path = base_path/'modeling'/'julia_stochastic_simulations'/'stochastic_sims.gzip'
    data = pd.DataFrame()
    if simulation_path.is_file(): 
        data = pd.read_parquet(simulation_path)

    # Rename columns/labels to match those for experimental data
    data.rename(columns={'copynum': 'copy_num', 'reg_gene': 'output', 'unreg_gene': 'marker'}, inplace=True)
    data['gene'] = data['design'].map({'Design 1': '1T', 'Design 2': '1T', 'Design 3': '1T',
                                            'Dual Vector': '2V', 'Dual Transcript': '2T'})
    data['design'] = data['design'].map({'Design 1': 1, 'Design 2': 2, 'Design 3': 3,
                                            'Dual Vector': 0, 'Dual Transcript': 0})
    data['kind'] = data['gene'] + '_' + data['design'].astype(str)
    data['condition'] = data['gene']+ '_' + data['moi'].astype(str)

    # Bin data and calculate statistics
    _, df_stats = calculate_bins_stats(data, by=['condition','design','moi','risc','gene','kind'])

    return data, df_stats

# later: move modeling results to different location
#   current base_path = rd.rootdir/'output'
def load_modeling_missplicing(base_path):

    # Load stochastic simulations
    simulation_path = base_path/'modeling'/'julia_missplicing_simulations'
    data_list = []
    for model in ['base', 'more_miRNA', 'more_mRNA']:
        if (simulation_path/('stochastic_sims_'+model+'.gzip')).is_file(): 
            d = pd.read_parquet(simulation_path/('stochastic_sims_'+model+'.gzip'))
            d['model'] = model
            data_list.append(d)
    
    # Combine data into a single dataframe
    data = pd.concat(data_list, ignore_index=True)

    # Rename columns/labels to match those for experimental data
    data.rename(columns={'copynum': 'copy_num', 'reg_gene': 'output', 'unreg_gene': 'marker'}, inplace=True)
    data['condition'] = data['model']+ '_' + data['moi'].astype(str)

    # Bin data and calculate statistics
    _, df_stats = calculate_bins_stats(data, by=['model','condition','risc','moi'])

    return data, df_stats

# later: move modeling results to different location
#   current base_path = rd.datadir/'projects'/'miR-iFFL'/'STRAIGHT-IN lines'
def load_ddpcr(base_path, metadata_path):

    # Load data
    ddpcr_path = rd.datadir/'projects'/'miR-iFFL'/'STRAIGHT-IN lines'
    ddpcr = pd.read_excel(ddpcr_path/'ddPCR.xlsx', header=1,)

    # Add metadata
    metadata = get_metadata(metadata_path/'construct-metadata.xlsx')
    ddpcr = ddpcr.merge(metadata, how='left', on='construct')

    return ddpcr
