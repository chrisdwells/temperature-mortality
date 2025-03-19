import pandas as pd
import numpy as np

"""
In Chen et al 2024 (https://www.nature.com/articles/s41467-024-45901-z), the
supplement states:

"Compared with the 0-64 age group, the 65-74 and ≥75 years age groups have 1.2 and 1.8 times
of heat-related mortality risks and 1.9 and 2.4 times of cold-related mortality risks""

We use that information to create scaling factors for modifying the total 
population mortality rate changes derived from the Bressler et al 2021 study.

Using 2015 age-stratified population fractions from 
https://population.un.org/dataportal/data/indicators/71/locations/900/start/
2015/end/2015/table/pivotbylocation?df=bc8978f2-c4e9-464f-8996-cbdd0154d9aa

M = total mortality
M_i = mortality of age group i
F_i = population fraction of group i

S_i,j = scaling mortality factor of i relative to j

M = sum(i) M_i*F_i (1)

M_i = M_j * S_i,j (2)

put (2) into (1)

M = sum(i) M_j*F_i*S_i,j = M_j * sum(i) F_i*S_i,j

so M_j/M = 1/(sum(i) * F_i*S_i,j) = 1/A_j (3)

with A_j = sum(i) * F_i*S_i,j

"""

# get the population fractions
pop_fractions = pd.read_csv('../data/unpopulation_dataportal_20250319111849.csv')
pop_fractions = pop_fractions[['Age', 'Value']]

pop_fractions_dict = {}

pop_fractions_dict['0-64'] = 0.01*(pop_fractions.loc[pop_fractions['Age'
                             ] == '0-14']['Value'].values[0
                         ] + pop_fractions.loc[pop_fractions['Age'
                             ] == '15-64']['Value'].values[0])

pop_fractions_dict['65-74'] = 0.01*(pop_fractions.loc[pop_fractions['Age'
                             ] == '65+']['Value'].values[0
                         ] - pop_fractions.loc[pop_fractions['Age'
                             ] == '75+']['Value'].values[0])

pop_fractions_dict['75+'] = 0.01*(pop_fractions.loc[pop_fractions['Age'
                             ] == '75+']['Value'].values[0
                         ])

# check population fractions sum to 1                                                         
assert np.around(pop_fractions_dict['0-64'] + pop_fractions_dict['65-74']
        + pop_fractions_dict['75+'], decimals=5) == 1                                                     


# organise the scaling factors
temps = ['Hot', 'Cold']
ages = ['0-64', '65-74', '75+']

# scaling_mortalities[temperature][top][bottom]
scaling_mortalities = {}
for temp in temps:
    scaling_mortalities[temp] = {}
    for age in ages:
        scaling_mortalities[temp][age] = {}

# from literature above
scaling_mortalities['Hot']['65-74']['0-64'] = 1.2
scaling_mortalities['Hot']['75+']['0-64'] = 1.8

scaling_mortalities['Cold']['65-74']['0-64'] = 1.9
scaling_mortalities['Cold']['75+']['0-64'] = 2.4

# create factors for each combination
for temp in temps:
    for a in ages:
        scaling_mortalities[temp][a][a] = 1
        for b in ages:
            if b not in scaling_mortalities[temp][a].keys():
                if a not in scaling_mortalities[temp][b].keys():
                    scaling_mortalities[temp][b][a] = scaling_mortalities[
                        temp][b]['0-64']/scaling_mortalities[
                            temp][a]['0-64']
                scaling_mortalities[temp][a][b] = 1/scaling_mortalities[temp][b][a]
                    
# calculate the coefficients as per equation (3)     
scaling_params = {}
for temp in temps:
    scaling_params[temp] = {}
    for j in ages:
        A_j = 0
        for i in ages:
            A_j += pop_fractions_dict[i]*scaling_mortalities[temp][i][j]
        scaling_params[temp][j] = 1/A_j
            
        #%%
# verify sum of population fractions * scaling factors is 1 (i.e so the 
# total mortality is unchanged upon age-stratifying)
for temp in temps:
    check_val = 0
    for age in ages:
        check_val += pop_fractions_dict[age]*scaling_params[temp][age]

    assert check_val == 1

#%%
# output to csv with FRIDA variables
# note names in FRIDA include the endpoint (e.g. 0-65 not 0-64) but
# they are consistent; should rename in FRIDA to be clear

frida_outputs = {} #pd.DataFrame(columns=['Variable', 'Value'])

frida_outputs['Demographics.Factor for hot mortality 0-65'
      ] = scaling_params['Hot']['0-64']
frida_outputs['Demographics.Factor for hot mortality 65-75'
      ] = scaling_params['Hot']['65-74']
frida_outputs['Demographics.Factor for hot mortality 75 plus'
      ] = scaling_params['Hot']['75+']

frida_outputs['Demographics.Factor for cold mortality 0-65'
      ] = scaling_params['Cold']['0-64']
frida_outputs['Demographics.Factor for cold mortality 65-75'
      ] = scaling_params['Cold']['65-74']
frida_outputs['Demographics.Factor for cold mortality 75 plus'
      ] = scaling_params['Cold']['75+']

df_outputs = pd.DataFrame(frida_outputs.items())

df_outputs.to_csv('../data/outputs/mortality_temperature_scalings_for_frida.csv', 
                     index=False)

