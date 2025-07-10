import numpy as np
import astropy.io.ascii
import pandas as pd

# all the targets in this analysis
targets = [
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
    "Per-emb-5",
    "Per-emb-12",
    "Per-emb-17",
    "Per-emb-18",
    "Per-emb-22",
    "Per-emb-27",
    "Per-emb-33",
    "Per-emb-35",
    "Per-emb-36",
    "Per-emb-44"
]

### PARSE `perseus1.txt` and `perseus2.txt` for PERSEUS

def ra_to_degrees(ra_str):
    h, m, s = map(float, ra_str.split(':'))
    return (h + m / 60 + s / 3600) * 15  # Convert hours to degrees

def dec_to_degrees(dec_str):
    sign = -1 if dec_str.startswith('-') else 1
    d, m, s = map(float, dec_str.lstrip('+-').split(':'))
    return sign * (d + m / 60 + s / 3600)

# read file
perseus1 = astropy.io.ascii.read("data/input/tobin2018_perseus_2.txt", delimiter="\t", guess=False).to_pandas()
# format name columns
perseus1['Main'] = perseus1['Source']
perseus1['Main'] = perseus1['Main'].apply(lambda x: str(x).removesuffix('-'+str(x).split('-')[-1]))
# filter for relevant targets
perseus1 = perseus1[perseus1['Main'].isin(targets)]
# order columns and rename
perseus1 = perseus1[['Main', 'Source', 'R.A.', 'Decl.']]
perseus1 = perseus1.rename(columns={'R.A.': 'RA', 'Decl.': 'Dec'})
# parse ra/dec columns
perseus1['RA'] = perseus1['RA'].apply(ra_to_degrees)
perseus1['Dec'] = perseus1['Dec'].apply(dec_to_degrees)
# set distance to 300pc
perseus1['Dis'] = 300

# read file
perseus2 = astropy.io.ascii.read("data/input/reynolds2024_perseus_1.txt", delimiter="\t", guess=False).to_pandas()
# format name columns
perseus2.loc[perseus2['Source'] == '-', 'Source'] = 'A'
perseus2['Source'] = perseus2['Source'].apply(lambda x: str(x).removeprefix('-'))
perseus2['Name'] = perseus2['Name'].replace('-', np.nan)
perseus2['Name'] = perseus2['Name'].fillna(method='ffill')
perseus2['Source'] = perseus2['Name'] + '-' + perseus2['Source']
# order columns and rename
perseus2 = perseus2[['Name', 'Source', 'RA', 'Dec']]
perseus2 = perseus2.rename(columns={'Name': 'Main'})
# filter for relevant targets
perseus2 = perseus2[perseus2['Main'].isin(targets)]
# parse ra/dec columns
perseus2['RA'] = perseus2['RA'].apply(ra_to_degrees)
perseus2['Dec'] = perseus2['Dec'].apply(dec_to_degrees)
# set distance to 300pc
perseus2['Dis'] = 300

# merge datasets
additional_targets = perseus2[perseus2['Main'].isin(['Per-emb-2', 'Per-emb-18'])]
perseus = pd.concat([perseus1, additional_targets]).sort_values('Source').reset_index(drop=True)
perseus['group'] = 'perseus'


### PARSE `distances.txt` for ORION

# read file
source_info = astropy.io.ascii.read("data/input/tobin2022_orion.txt").to_pandas()

# filter for relevant targets
source_info = source_info[source_info['Main'].isin(targets)]

# create ra/dec columns
source_info['RA'] = 15 * (source_info['RAh'] + source_info['RAm'] / 60 + source_info['RAs'] / 3600)
source_info['DE-'] = source_info["DE-"].apply(lambda x: -1 if x == "-" else 1)
source_info['Dec'] = source_info['DE-'] * (source_info['DEd'] + source_info['DEm'] / 60 + source_info['DEs'] / 3600)

# order columns and sort
source_info = source_info[['Main', 'Source', 'RA', 'Dec', 'Dis', 'LBol', 'TBol', 'Class', 'SigmaYSO']]
source_info = source_info.sort_values(['Main', 'Source']).reset_index(drop=True)

# fix HOPS-361-N
source_info["Main"] = source_info["Main"].replace({"HOPS-361": "HOPS-361-N"})
source_info["Source"] = source_info["Source"].apply(lambda x: str(x).replace("HOPS-361", "HOPS-361-N"))

source_info['group'] = 'orion'

# save
source_info = pd.concat([source_info, perseus])
source_info.to_csv("data/output/source_info.csv",index=False)



### PARSE `notes.txt`

# used to parse the `Binary` column
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

# used to parse the `outflow_source` column
def fix_outflow_source(x):
    if str(x).__contains__('+'):
        return 'both'
    elif len(str(x)) > 1:
        return str(x)[0] + '-' + str(x)[1]
    else:
        return x
    
def angle_west_of_north(v):
    x, y = v
    angle_rad = np.arctan2(y, x)  # atan2(y, x) gives the angle in radians
    angle_deg = np.degrees(angle_rad)  # Convert to degrees
    return (90 - angle_deg) % 360  # Ensure the angle is in [0, 360) range

# read file
df = pd.read_csv("data/input/notes.csv")
# keep relevant columns
df = df[['Field', 'Binary', 'Outflow Source', 'Blue Channels', 'Red Channels', 'Red Center Corrected', 'Blue Center Corrected', 'Average Angle (Red)', 'Lobe PA Offset']]
# get sources with secondaries
no_binary = df['Binary'].isna()
# create columns to identify each source in the pair
df.loc[~no_binary, 'source_a'] = df.loc[~no_binary, 'Binary'].apply(lambda x: split_source_string(x)[0])
df.loc[~no_binary, 'source_b'] = df.loc[~no_binary, 'Binary'].apply(lambda x: split_source_string(x)[1])
df.loc[no_binary, 'source_a'] = df.loc[no_binary, 'Outflow Source']
df.loc[no_binary, 'source_b'] = None
# rename columns
df = df[['Field', 'source_a', 'source_b', 'Outflow Source', 'Blue Channels', 'Red Channels', 'Blue Center Corrected', 'Red Center Corrected', 'Average Angle (Red)', 'Lobe PA Offset']].rename(columns={'Field': 'field', 'Average Angle (Red)': 'outflow_PA', 'Outflow Source': 'outflow_source', 'Red Channels': 'red_channels', 'Blue Channels': 'blue_channels', 'Red Center Corrected': 'red_outflow_PA', 'Blue Center Corrected': 'blue_outflow_PA', 'Lobe PA Offset': 'lobe_PA_offset'})

# merge coordinates and distance
new_rows = []
for i, row in df.iterrows():

    ### exceptions
    # removal of Per-emb-5, not a binary
    if row['field'] == 'Per-emb-5':
        continue
    # # removal of secondary outflows on stars >500 AU from inner binary
    # if (row['field'] == 'HOPS-12') & (pd.isna(row['source_b'])):
    #     continue
    # if (row['field'] == 'HOPS-203') & (pd.isna(row['source_b'])):
    #     continue

    key_a = row['field'] + '-' + row['source_a']
    row_a = source_info.loc[source_info['Source'] == key_a]
    row['source_a_ra'] = row_a['RA'].iloc[0]
    row['source_a_dec'] = row_a['Dec'].iloc[0]
    row['group'] = row_a['group'].iloc[0]

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
master = pd.DataFrame(new_rows)[['group', 'field', 'source_a', 'source_a_ra', 'source_a_dec', 'source_b', 'source_b_ra', 'source_b_dec', 'distance', 'outflow_source', 'red_channels', 'blue_channels', 'red_outflow_PA', 'blue_outflow_PA', 'outflow_PA', 'lobe_PA_offset', 'binary_PA']]
# parse `outflow_source`
master['outflow_source'] = master['outflow_source'].apply(lambda x: fix_outflow_source(x))

# compute delta_PA
angle = np.abs(master['outflow_PA'] - master['binary_PA']) % 180
angle = np.min([angle, 180 - angle], axis=0)
master['delta_PA'] = angle

# save
master.to_csv('data/output/outflow_data.csv',index=False)


### PARSE Separations Data

ori = astropy.io.ascii.read("data/input/tobin2022_orion_pairings.txt").to_pandas()
per = astropy.io.ascii.read("data/input/tobin2022_perseus_pairings.txt").to_pandas()

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
}

sep = pd.concat([ori[ori['Pair'].isin(list(ori_keys.keys()))], per[per['Pair'].isin(list(per_keys.keys()))]]).reset_index(drop=True)
sep['Pair'] = sep['Pair'].apply(lambda x: {**per_keys, **ori_keys}[x])
sep['Tertiary'] = sep['Class'].str.contains('\\(')

def class_map(x):
    if 'C0' in x:
        return 'C0'
    elif 'CI' in x:
        return 'C1'
    elif 'FS' in x:
        return 'FS'
    else:
        return np.nan

sep['Class'] = sep['Class'].apply(class_map)
sep = sep[['Pair', 'Tertiary', 'PSep', 'Class']].sort_values(by='Pair').reset_index(drop=True)
sep = sep.rename(columns={'PSep': 'separation', 'Class': 'class', 'Tertiary': 'tertiary', 'Pair': 'pair'})
sep.to_csv("data/output/separations.csv",index=False)