"""
Description
-----------

Convert density-dependent experiment data to a format in which we can estimate
the transmission and the growth function.

Author: Mark Wilber
"""

import numpy as np
import pandas as pd

dat = pd.read_csv("../archival/DensExpt.csv")
dat.set_index("swab_id", inplace=True)
dat.ix["MDF16_13_Jul31", "ZE"] = 0 # Replacing missing data 
dat.reset_index(inplace=True)

# If day is absent drop that row from the data
dat = dat[~np.isnan(dat.Day)]

# Add some extra columns to data for easier subsetting later
tank_info = dat.groupby("Tank_num").agg({'InitialFrogs' :
                                                lambda x: np.unique(x)[0]})
dat_new = dat.set_index("Tank_num").join(tank_info, rsuffix="_full").reset_index()

# Specify a swab_event to classify tadpoles and frogs by swab event
class_dict = {1: [0, -3], 2: range(3, 6), 3: range(9, 12), 4: range(16, 21),
              5: range(23, 26), 6: range(30, 33), 7: range(36, 39),
              8: range(45, 48), 9: range(52, 55), 10: range(60, 63),
              11: range(66, 69), 12: range(73, 76), 13: range(136, 139)}

# Make a swab dictionary for varying dates
swab_dict = {}
for key, value in class_dict.iteritems():
    for v in value:
        swab_dict[int(v)] = key

# Make a column that describes the swab event
dat_new['swab_event'] = [swab_dict[int(day)] for day in dat_new.Day]
dat_new['group_id'] = ["_".join([str(int(a)), str(int(b)), str(int(c))]) for a, b, c in
                    np.array(dat_new[['Tank_num', 'InitialFrogs_full', 'swab_event']])]

# Subset data on frog and tadpoles
ind = dat_new.Frog_or_Tadpole == 'F'
frog_dat = dat_new[ind]
tad_dat = dat_new[~ind]

# Find the sum of zoospores for the tadpoles in each tank at a given time point
sum_zes = tad_dat.groupby('group_id').agg({'ZE': np.sum})
frog_dat = frog_dat.set_index("group_id").join(sum_zes, rsuffix="_tadpoles").reset_index()

# Find the number of infected Frogs and their total zoospore load at a given
# time point
sum_frogs_zes = frog_dat.groupby("group_id").agg({'ZE': np.nansum})
num_frogs_inf = frog_dat.groupby("group_id").agg({'ZE': lambda x: np.nansum(x != 0)})

inf_info = num_frogs_inf.rename(columns={'ZE': "num_infected"})
inf_info['ze_frogs'] = sum_frogs_zes.ZE

# Include the frog infectedness and ZE load into the full matrix
frog_dat = frog_dat.set_index("group_id").join(inf_info).reset_index()

## Reformat the data for analysis ##

# Split by individual
individuals = frog_dat.Indiv_id.unique()

ind_dfs = []

for indiv in individuals:

    tdat = frog_dat[frog_dat.Indiv_id == indiv]

    # Sort by Day
    tdat = tdat.sort_values(by="Day")

    time_diff = np.array(tdat.Day[1:]) - np.array(tdat.Day[0:-1])

    matched_ZE = zip(tdat.ZE[0:-1], tdat.ZE[1:], tdat.swab_id[1:],
                            time_diff, tdat.Day[1:],
                            np.repeat(indiv, len(time_diff)),
                            tdat.InitialFrogs[:-1],
                            tdat.weight[:-1], tdat.svl[:-1],
                            tdat.sex[:-1], tdat.ZE_tadpoles[:-1],
                            tdat.Tank_num[:-1], tdat.ze_frogs[:-1],
                            tdat.num_infected[:-1])

    tdf = pd.DataFrame(matched_ZE, columns=['size', 'sizeNext', 'swab_id', 'time',
                                "Day", 'individual', 'density',
                                "weight", 'svl', 'sex', 'ze_tadpoles',
                                'tank_num', "ze_frogs", "num_infected"])

    tdf['cum_ze_frogs'] = np.cumsum(tdf.ze_frogs)
    tdf['cum_ze_tadpoles'] = np.cumsum(tdf.ze_tadpoles)

    ind_dfs.append(tdf)

full_data = pd.concat(ind_dfs)

# Drop the points where they are dead in both time points
full_data = full_data[~np.isnan(full_data['size'])]

# Save the converted data
full_data.to_csv("formatted_dd_data.csv", index=False)

