import numpy as np
import astropy.io.ascii
import pandas as pd
from astropy.table import Table
from io import StringIO
from astropy import units as u
from _create_figs import getIdx
import os
import io
import glob
import yaml
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.io.fits as fits
import astropy.constants as const

### SCRIPT OPTIONS
# output
output_folder = "results/tables"
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
if not os.path.exists("data/output"):
    os.mkdir("data/output")

### SCRIPT
# all the targets in this analysis
ALL_TARGETS = [
    "HH270VLA1",
    "HOPS-323",
    "HOPS-312",
    "HOPS-364",
    "HOPS-45",
    "HOPS-361", # file is named 361-N
    "HOPS-366", 
    "HOPS-384",
    "HOPS-304",
    "HOPS-357",
    "HOPS-363",
    "HOPS-400",
    "HOPS-290",
    "HOPS-56",
    "HOPS-193",
    "HOPS-242",
    "HOPS-70",
    "HOPS-138",
    "HOPS-85",
    "HOPS-182",
    "HOPS-28",
    "HOPS-92",
    "HOPS-32",
    "HOPS-158",
    "HOPS-75",
    "HOPS-395",
    "HOPS-43", # where did this one come from?
    "HOPS-255",
    "HOPS-213",
    "HOPS-77",
    "HOPS-281",
    "HOPS-173",
    "HOPS-282",
    "HOPS-84",
    "HOPS-288",
    "HOPS-248",
    "HOPS-168",
    "HOPS-203",
    "HOPS-163",
    "HOPS-12",
    "HOPS-261",

    "Per-emb-2",
    "Per-emb-12",
    "Per-emb-17",
    "Per-emb-18",
    "Per-emb-22",
    "Per-emb-27",
    "Per-emb-33",
    "Per-emb-35",
    "Per-emb-36",
    "Per-emb-44",
    "L1448 IRS3C"
]

def write_latex_table(df, filename):
    table = Table.from_pandas(df)
    latex_buffer = io.StringIO()
    table.write(
        latex_buffer,
        format='ascii.latex',
        latexdict={'tabletype': 'deluxetable'},
        overwrite=True
    )
    latex_string = latex_buffer.getvalue().replace('  ', r' \nodata ')
    with open(f'results/tables/{filename}', 'w') as f:
        f.write(latex_string)

# used to parse perseus RA/Dec columns
def ra_to_degrees(ra_str):
    h, m, s = map(float, ra_str.split(':'))
    return (h + m / 60 + s / 3600) * 15  # Convert hours to degrees
def dec_to_degrees(dec_str):
    sign = -1 if dec_str.startswith('-') else 1
    d, m, s = map(float, dec_str.lstrip('+-').split(':'))
    return sign * (d + m / 60 + s / 3600)

# used to parse the `Binary` in `notes.csv`
def split_source_string(x):
    x = str(x).split('+')
    if len(x[0]) == 1:
        a = x[0]
    else:
        a = x[0][0] + '-' + x[0][1]
    if len(x[1]) == 1:
        b = x[1]
    else:
        b = x[1][0] + '-' + x[1][1]
    return [a, b]

# used to parse the `outflow_source` column in `notes.csv`
def fix_outflow_source(x):
    if str(x).__contains__('+'):
        return 'both'
    elif len(str(x)) > 1:
        return str(x)[0] + '-' + str(x)[1]
    else:
        return x
    
def angle_west_of_north(v):
    x, y = v
    angle_rad = np.arctan2(y, x)
    angle_deg = np.degrees(angle_rad)
    return (90 - angle_deg) % 360 

# used to rename the `Class` column in `separations.csv`
def class_map(x):
    if 'C0' in x:
        return '0'
    elif 'CI' in x:
        return 'I'
    elif 'FS' in x:
        return 'FS'
    else:
        return np.nan


### PARSE `tobin2018_perseus_2.txt` and `reynolds2024_perseus_1.txt` for PERSEUS
# read file
tobin2018_perseus_2 = astropy.io.ascii.read("data/input/tobin2018_perseus_2.txt", delimiter="\t", guess=False).to_pandas()
# format name columns
tobin2018_perseus_2['Main'] = tobin2018_perseus_2['Source'].apply(lambda x: str(x).removesuffix('-'+str(x).split('-')[-1]))
# select columns and rename
tobin2018_perseus_2 = tobin2018_perseus_2[['Main', 'Source', 'R.A.', 'Decl.']].rename(columns={'R.A.': 'RA', 'Decl.': 'Dec'})
# read file
reynolds2024_perseus_1 = astropy.io.ascii.read("data/input/reynolds2024_perseus_1.txt", delimiter="\t", guess=False).to_pandas()
# format name columns
reynolds2024_perseus_1.loc[reynolds2024_perseus_1['Source'] == '-', 'Source'] = 'A'
reynolds2024_perseus_1['Name'] = reynolds2024_perseus_1['Name'].replace('-', np.nan).fillna(method='ffill')
reynolds2024_perseus_1['Source'] = reynolds2024_perseus_1['Source'].apply(lambda x: str(x).removeprefix('-'))
reynolds2024_perseus_1['Source'] = reynolds2024_perseus_1['Name'] + '-' + reynolds2024_perseus_1['Source']
# select columns and rename
reynolds2024_perseus_1 = reynolds2024_perseus_1[['Name', 'Source', 'RA', 'Dec', 'T _B']].rename(columns={'Name': 'Main', 'T _B': 'TBol'})

### MERGE PERSEUS datasets into `perseus`
perseus = pd.merge(tobin2018_perseus_2, reynolds2024_perseus_1, on=['Main', 'Source'], how='outer', suffixes=(None, '_y'))
# filter for relevant targets
perseus = perseus[perseus['Main'].isin(ALL_TARGETS)]
# parse ra/dec columns
perseus['RA'] = perseus['RA'].fillna(perseus['RA_y'])
perseus['Dec'] = perseus['Dec'].fillna(perseus['Dec_y'])
perseus = perseus.drop(columns=['RA_y', 'Dec_y'])
perseus['RA'] = perseus['RA'].apply(ra_to_degrees)
perseus['Dec'] = perseus['Dec'].apply(dec_to_degrees)
# set distance to 300pc
perseus['Dis'] = 300
# add group column
perseus['group'] = 'perseus'

### PARSE `distances.txt` for ORION
# read file
orion = astropy.io.ascii.read("data/input/tobin2022_orion.txt").to_pandas()
# filter for relevant targets
orion = orion[orion['Main'].isin(ALL_TARGETS)]
# create ra/dec columns
orion['RA'] = 15 * (orion['RAh'] + orion['RAm'] / 60 + orion['RAs'] / 3600)
orion['DE-'] = orion["DE-"].apply(lambda x: -1 if x == "-" else 1)
orion['Dec'] = orion['DE-'] * (orion['DEd'] + orion['DEm'] / 60 + orion['DEs'] / 3600)
# select columns and sort
orion = orion[['Main', 'Source', 'RA', 'Dec', 'Dis', 'LBol', 'TBol', 'Class', 'SigmaYSO']]
# fix HOPS-361-N
orion["Main"] = orion["Main"].replace({"HOPS-361": "HOPS-361-N"})
orion["Source"] = orion["Source"].apply(lambda x: str(x).replace("HOPS-361", "HOPS-361-N"))
# add group column
orion['group'] = 'orion'

### MERGE ORION and PERSEUS into `source_info`
source_info = pd.concat([orion, perseus]).sort_values(by=['RA']).reset_index(drop=True)
source_info.to_csv("data/output/source_info.csv",index=False)

### PARSE `notes.txt`
# read file
df = pd.read_csv("data/input/notes.csv")
# get sources with secondaries
no_binary = df['Binary'].isna()
# create columns to identify each source in the pair
df.loc[~no_binary, 'source_a'] = df.loc[~no_binary, 'Binary'].apply(lambda x: split_source_string(x)[0])
df.loc[~no_binary, 'source_b'] = df.loc[~no_binary, 'Binary'].apply(lambda x: split_source_string(x)[1])
df.loc[no_binary, 'source_a'] = df.loc[no_binary, 'Outflow Source']
df.loc[no_binary, 'source_b'] = None
# select and rename columns
df = df[['Field', 'Confidence', 'source_a', 'source_b', 'Outflow Source', 'Blue Channels', 'Red Channels', 'Blue Center Corrected', 'Red Center Corrected', 'Average Angle (Blue)', 'Lobe PA Offset']].rename(columns={'Field': 'field', 'Average Angle (Blue)': 'outflow_PA', 'Outflow Source': 'outflow_source', 'Red Channels': 'red_channels', 'Blue Channels': 'blue_channels', 'Red Center Corrected': 'red_outflow_PA', 'Blue Center Corrected': 'blue_outflow_PA', 'Lobe PA Offset': 'lobe_PA_offset', 'Confidence': 'confidence'})

# merge coordinates and distance
new_rows = []
for i, row in df.iterrows():
    # construct source_a
    key_a = row['field'] + '-' + row['source_a']
    row_a = source_info.loc[source_info['Source'] == key_a]
    row['source_a_ra'] = row_a['RA'].iloc[0]
    row['source_a_dec'] = row_a['Dec'].iloc[0]
    row['group'] = row_a['group'].iloc[0]

    # construct source_b, if source is not a tertiary
    if row['source_b'] != None:
        key_b = row['field'] + '-' + row['source_b']
        row_b = source_info.loc[source_info['Source'] == key_b]
        row['source_b_ra'] = row_b['RA'].iloc[0]
        row['source_b_dec'] = row_b['Dec'].iloc[0]

        # calculate angle of the separation vector
        separation_vector = np.array([row['source_b_ra'] - row['source_a_ra'], row['source_b_dec'] - row['source_a_dec']])
        row['binary_PA'] = angle_west_of_north(separation_vector)
    else:
        row['source_b_ra'] = None
        row['source_b_dec'] = None
        row['binary_PA'] = None

    row['distance'] = row_a['Dis'].iloc[0]
    new_rows.append(row)

# create new dataframe
outflow_data = pd.DataFrame(new_rows)[['group', 'field', 'confidence', 'source_a', 'source_a_ra', 'source_a_dec', 'source_b', 'source_b_ra', 'source_b_dec', 'distance', 'outflow_source', 'red_channels', 'blue_channels', 'red_outflow_PA', 'blue_outflow_PA', 'outflow_PA', 'lobe_PA_offset', 'binary_PA']]
# parse `outflow_source`
outflow_data['outflow_source'] = outflow_data['outflow_source'].apply(lambda x: fix_outflow_source(x))
# compute delta_PA
angle = np.abs(outflow_data['outflow_PA'] - outflow_data['binary_PA']) % 180
angle = np.min([angle, 180 - angle], axis=0)
outflow_data['delta_PA'] = angle
# sort
outflow_data = outflow_data.sort_values('source_a_ra')
# save
outflow_data.to_csv('data/output/outflow_data.csv',index=False)


### PARSE Separations Data
# read files
orion_sep = astropy.io.ascii.read("data/input/tobin2022_orion_pairings.txt").to_pandas()
perseus_sep = astropy.io.ascii.read("data/input/tobin2022_perseus_pairings.txt").to_pandas()
ori_keys = {
    'H12-B-A+H12-B-B' : 'HOPS-12', # double
    '(H12-B-A+H12-B-B)+H12-A' : 'HOPS-12',
    'H290-A+H290-B' : 'HOPS-290',
    'H92-A-B+H92-A-A' : 'HOPS-92', # double
    'H92-B+(H92-A-B+H92-A-A)' : 'HOPS-92',
    'H288-A-B+H288-A-A' : 'HOPS-288', # double
    'H288-B+(H288-A-B+H288-A-A)' : 'HOPS-288',
    'H203-B+H203-A' : 'HOPS-203', # double
    'H203-C+(H203-B+H203-A)' : 'HOPS-203',
    'HH270VLA1-B+HH270VLA1-A': 'HH270VLA1',
    'H32-B+H32-A' : 'HOPS-32',
    'H84-B+H84-A' : 'HOPS-84',
    'H168-A+H168-B' : 'HOPS-168',
    'H281-B+H281-A' : 'HOPS-281',
    'H312-B+H312-A' : 'HOPS-312',
    'H364-A+H364-B' : 'HOPS-364',
    'H395-B+H395-A' : 'HOPS-395',
    'H400-B+H400-A' : 'HOPS-400',
    'H182-B+H182-A' : 'HOPS-182',
    'H323-B+H323-A' : 'HOPS-323',
    'H282-B+H282-A' : 'HOPS-282',
    'H366-A+H366-B' : 'HOPS-366',
    'H193-B+H193-A' : 'HOPS-193',
    'H304-B+H304-A' : 'HOPS-304',
    'H361-C-A+H361-C-B' : 'HOPS-361-N',
    'H384-A+H384-A-B' : 'HOPS-384',
    'H363-B+H363-A' : 'HOPS-363',
    'H173-B+H173-A' : 'HOPS-173',
    'H75-B+H75-A' : 'HOPS-75',
    'H213-A+H213-B' : 'HOPS-213',
    'H158-B+H158-A' : 'HOPS-158'
}
per_keys = {
    'P2-A-P2-B' : 'Per-emb-2',
    'P12-A-P12-B' : 'Per-emb-12',
    'P17-A-P17-B' : 'Per-emb-17',
    'P18-A-P18-B' : 'Per-emb-18',
    'P22-A-P22-B' : 'Per-emb-22',
    'P27-A-P27-B' : 'Per-emb-27',
    'P33-B-P33-C' : 'Per-emb-33', # double
    'P33-A-(P33-B+P33-C)' : 'Per-emb-33',
    'P35-A-P35-B' : 'Per-emb-35',
    'P36-A-P36-B' : 'Per-emb-36',
    'P44-A-P44-B' : 'Per-emb-44',
    'L1448NW-A-L1448NW-B' : 'L1448 IRS3C',
}
# merge datasets and filter for relevant targets
sep = pd.concat([orion_sep[orion_sep['Pair'].isin(list(ori_keys.keys()))], perseus_sep[perseus_sep['Pair'].isin(list(per_keys.keys()))]]).reset_index(drop=True)
# rename pairs
sep['Pair'] = sep['Pair'].apply(lambda x: {**per_keys, **ori_keys}[x])
# note tertiary sources
sep['Tertiary'] = sep['Class'].str.contains('\\(')

sep['Class'] = sep['Class'].apply(class_map)
sep = sep[['Pair', 'Tertiary', 'PSep', 'Class']].rename(columns={'PSep': 'separation', 'Class': 'class', 'Tertiary': 'tertiary', 'Pair': 'pair'})
sep = sep.sort_values(by='pair').reset_index(drop=True)
# sep.to_csv("data/output/separations.csv",index=False)

# CREATE temporary TABLE with velocity ranges and field coordinates from HDU DATA
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)
IMAGE_DIRECTORY = config["data_dir"]
files = glob.glob(f'{IMAGE_DIRECTORY}/*/*12co*.fits') + glob.glob(f'{IMAGE_DIRECTORY}/*/*spw39*.fits')

def channel_to_velocity(wcs, channel):
    pix = np.array([[0,0,channel]])
    world = wcs.all_pix2world(pix, 0)
    _, _, freq = world[0]
    rest_freq = 230.538 * u.GHz
    freq = freq * u.Hz
    v = (rest_freq - freq) / rest_freq * const.c
    return round(v.to(u.km / u.s).value, 4)

hdu_data = []
for file in files:
    target_name = os.path.dirname(file).split('/')[-1]
    hdu = fits.open(file)[0]
    wcs = WCS(hdu.header, naxis=3)
    center = SkyCoord(hdu.header['OBSRA'], hdu.header['OBSDEC'], unit=u.degree)

    if target_name in outflow_data['field'].apply(str.casefold).values:
        for i, row in outflow_data.loc[outflow_data['field'].apply(str.casefold) == target_name].iterrows():
            red_channels = getIdx([row['red_channels']])
            blue_channels = getIdx([row['blue_channels']])

            if len(red_channels) != 0:
                red_v_start = channel_to_velocity(wcs, red_channels[0])
                red_v_end = channel_to_velocity(wcs, red_channels[-1])
                red_v_range = f"{red_v_start} to {red_v_end}"
            else:
                red_v_range = np.nan
            

            if len(blue_channels) != 0:
                blue_v_start = channel_to_velocity(wcs, blue_channels[0])
                blue_v_end = channel_to_velocity(wcs, blue_channels[-1])
                blue_v_range = f"{blue_v_start} to {blue_v_end}"
            else:
                blue_v_range = np.nan
            
            d = {'field2': target_name, 'source_a_ra': row['source_a_ra'], 'ra': center.ra.value, 'dec': center.dec.value, 'red_v_range': red_v_range, 'blue_v_range': blue_v_range}
            hdu_data.append(d)
    else:
        d = {'field2': target_name, 'source_a_ra': np.nan, 'ra': center.ra.value, 'dec': center.dec.value, 'red_v_range': np.nan, 'blue_v_range': np.nan}
        hdu_data.append(d)

df_hdudata = pd.DataFrame(hdu_data).sort_values('source_a_ra').reset_index(drop=True)
df_hdudata = df_hdudata.groupby(['field2']).agg('first').reset_index()

### CREATE by_outflow TABLE
# note to self this is not creating a copy, just an alias
by_outflow = outflow_data.copy()
# combine sources for 'both' outflow sources
by_outflow.loc[by_outflow['outflow_source'] == 'both', 'outflow_source'] = by_outflow.loc[by_outflow['outflow_source'] == 'both']['source_a'] + '+' + by_outflow.loc[by_outflow['outflow_source'] == 'both']['source_b']
# note tertiary sources
by_outflow['tertiary'] = by_outflow['binary_PA'].isna()
# merge separations
by_outflow = pd.merge(by_outflow, sep, left_on=['field', 'tertiary'], right_on=['pair', 'tertiary'], how='left')
# define type column
by_outflow['type'] = by_outflow['tertiary'].map({True: 'Tertiary', False: 'Binary'})
# select columns
by_outflow = by_outflow[['field', 'confidence', 'source_a_ra', 'type', 'separation', 'class', 'outflow_source', 'outflow_PA', 'lobe_PA_offset', 'binary_PA', 'delta_PA']]
# format missing values
by_outflow = by_outflow.fillna('-')
by_outflow[['outflow_PA', 'lobe_PA_offset', 'binary_PA', 'delta_PA']] = by_outflow[['outflow_PA', 'lobe_PA_offset', 'binary_PA', 'delta_PA']].replace('-', np.nan).astype(float).round(0).astype('Int64')
# sort dataframe
by_outflow = by_outflow.sort_values('source_a_ra').drop('source_a_ra',axis=1).reset_index(drop=True)
# save to CSV
by_outflow.to_csv('data/output/data_by_outflow.csv', index=False)
# save to LaTeX table
write_latex_table(by_outflow, 'by-outflow.tex')

### CREATE by_field TABLE
# get info for sources without outflows
unmeasured_sources = ['HOPS-56', 'HOPS-70', 'HOPS-261', 'HOPS-158', 'HOPS-138', 'HOPS-45', 'HOPS-85', 'HOPS-77'] + ['HOPS-28', 'HOPS-163', 'HOPS-242', 'HOPS-248', 'HOPS-255', 'HOPS-357']
sources_no_outflows = source_info[['Main', 'Source', 'RA', 'Dec']].groupby('Main').agg({
    'Source': 'count',
    'RA': 'first',
    'Dec': 'first'
}).reset_index().rename(columns={'Source': 'n_sources', 'RA': 'ra', 'Dec': 'dec'})
# extract number of sources for all fields for later use
n_sources = sources_no_outflows[['Main', 'n_sources']]
# filter for sources without outflows
sources_no_outflows = sources_no_outflows.loc[sources_no_outflows['Main'].isin(unmeasured_sources)].drop(columns='n_sources').reset_index(drop=True)

# define other names for fields
other_names = pd.DataFrame({
    "Per-emb-2": "IRAS 03292+3039",
    "Per-emb-12": "NGC 1333 IRAS4A",
    "Per-emb-22": "L1448 IRS2",
    "Per-emb-27": "NGC 1333 IRS2A",
    "Per-emb-33": "L1448 IRS3B",
    "Per-emb-36": "NGC 1333 IRAS2B",
    "Per-emb-44": "SVS 13A",
    "L1448 IRS3C": "Per-emb-107",
}, index=[0]).T.reset_index().rename(columns={'index': 'field', 0: 'other_names'})

# start by grouping outflow data by field
by_field = outflow_data.groupby(['field']).agg('first').reset_index()
# select and rename columns
by_field = by_field[['field', 'source_a_ra', 'source_a_dec', 'red_channels', 'blue_channels']].rename(columns={'source_a_ra': 'ra', 'source_a_dec': 'dec'})
# add in unmeasured sources
by_field = pd.concat([by_field, sources_no_outflows.rename(columns={'Main': 'field'})], ignore_index=True)
# add in n_sources column
by_field = by_field.merge(n_sources, left_on='field', right_on='Main', how='left')
# add other names column
by_field = by_field.merge(other_names, on='field', how='left')
# merge ra, dec, red_v, blue_v from hdudata
by_field['field2'] = by_field['field'].apply(str.casefold)
by_field = by_field.merge(df_hdudata, on='field2', suffixes=['_x', None])[['field', 'other_names', 'n_sources', 'ra', 'dec', 'red_v_range', 'blue_v_range']]
# combine channel columns
by_field['integrated_channels'] = (by_field['red_v_range'].fillna('') + ', ' + by_field['blue_v_range'].fillna('')).str.lstrip(', ').str.rstrip(', ')
# remove duplicates
by_field = by_field.drop_duplicates(subset=['field']).reset_index(drop=True)
# format number columns
by_field[['ra', 'dec']] = by_field[['ra', 'dec']].round(5)
by_field[['n_sources']] = by_field[['n_sources']].round(0).astype(int)
# order columns and sort
by_field = by_field[['field', 'other_names', 'n_sources', 'ra', 'dec', 'integrated_channels']].sort_values('ra').reset_index(drop=True)
# save to CSV
by_field.to_csv('data/output/data_by_field.csv', index=False)
# save to LaTeX table
write_latex_table(by_field, 'by-field.tex')