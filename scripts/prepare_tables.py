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
    "HOPS-43",
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
    "HOPS-261"
]

### PARSE `distances.txt`

# read file
source_info = astropy.io.ascii.read("data/input/distances.txt").to_pandas()

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

# save
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
df = df[['Field', 'Binary', 'Outflow Source', 'Blue Channels', 'Red Channels', 'Average Angle (Blue)']]
# remove sources with secondaries (deal with those later)
df = df[~df['Binary'].isna()]
# create columns to identify each source in the pair
df['source_a'] = df['Binary'].apply(lambda x: split_source_string(x)[0])
df['source_b'] = df['Binary'].apply(lambda x: split_source_string(x)[1])
# rename columns
df = df[['Field', 'source_a', 'source_b', 'Outflow Source', 'Blue Channels', 'Red Channels', 'Average Angle (Blue)']].rename(columns={'Field': 'field', 'Average Angle (Blue)': 'outflow_angle', 'Outflow Source': 'outflow_source', 'Red Channels': 'red_channels', 'Blue Channels': 'blue_channels'})

# merge coordinates and distance
new_rows = []
for i, row in df.iterrows():
    key_a = row['field'] + '-' + row['source_a']
    key_b = row['field'] + '-' + row['source_b']
    row_a = source_info.loc[source_info['Source'] == key_a]
    row_b = source_info.loc[source_info['Source'] == key_b]
    row['source_a_ra'] = row_a['RA'].iloc[0]
    row['source_a_dec'] = row_a['Dec'].iloc[0]
    row['source_b_ra'] = row_b['RA'].iloc[0]
    row['source_b_dec'] = row_b['Dec'].iloc[0]
    row['distance'] = row_a['Dis'].iloc[0]
    separation_vector = np.array([row['source_b_ra'] - row['source_a_ra'], row['source_b_dec'] - row['source_a_dec']])
    row['separation_angle'] = angle_west_of_north(separation_vector)
    new_rows.append(row)

# create new dataframe
master = pd.DataFrame(new_rows)[['field', 'source_a', 'source_a_ra', 'source_a_dec', 'source_b', 'source_b_ra', 'source_b_dec', 'distance', 'outflow_source', 'red_channels', 'blue_channels', 'outflow_angle', 'separation_angle']]
# parse `outflow_source`
master['outflow_source'] = master['outflow_source'].apply(lambda x: fix_outflow_source(x))
# save
master.to_csv('data/output/outflow_data.csv',index=False)