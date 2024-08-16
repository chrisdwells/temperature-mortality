import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import copy
import scipy.stats

def fit(x, a, beta_t, beta_t2):
    yhat = a + beta_t*x + beta_t2*x**2
    return yhat

impacts = {
    'hot':'red', 
    'cold':'blue', 
    'net':'black',
    }

percs = ['low', 'mid', 'high']

gmst = 'data/AR6_GMST.csv'
gmst_in = pd.read_csv(gmst)

gmst_2000_2019 = np.mean(gmst_in[(2000 <= gmst_in['year']) & 
                     (gmst_in['year']<= 2019)]['gmst'])

gmst_offset = gmst_2000_2019
# gmst_1970_1999 = np.mean(gmst_in[(1970 <= gmst_in['year']) & 
#                      (gmst_in['year']<= 1999)]['gmst'])

# gmst_offset = gmst_2000_2019 - gmst_1970_1999

#%%

temps_plot = np.linspace(0, 4.5, 100)

params = {}

for impact in impacts.keys():
    params[impact] = {}
    
    for perc in percs:
            
        mort_file = f'data/Bressler_fig3_{impact}_{perc}.csv'
        df_in = pd.read_csv(mort_file)
        
        temps = df_in['Temperature'].values
        morts = df_in['Mortality'].values
    
        temps_offset = temps + gmst_offset
    
        plt.scatter(temps_offset, morts, color=impacts[impact])
    
    
        params_in, _ = curve_fit(
            fit, temps_offset, morts)
        
        params_no_int = copy.deepcopy(params_in)
    
        params_no_int[0] = 0
        
        params[impact][perc] = params_no_int
    
        if perc == 'mid':
                
            plt.plot(temps_plot, fit(temps_plot, *params_no_int), color=impacts[impact],
                     label=f'{impact}')
        else:
            
            plt.plot(temps_plot, fit(temps_plot, *params_no_int), color=impacts[impact])
        
        plt.plot(temps_plot, fit(temps_plot, *params_in), color=impacts[impact],
                 linestyle='dashed')

plt.xlabel(f'GMST cf pi (offset: {np.around(gmst_offset, decimals=3)}K)')
plt.ylabel('Mortality cf pi (solid), cf 2010 (dashed)')
plt.legend()

#%%

def opt(x, q025_desired, q50_desired, q975_desired):
    q025, q50, q975 = scipy.stats.skewnorm.ppf(
        (0.025, 0.50, 0.975), x[0], loc=x[1], scale=x[2]
    )
    return (q025 - q025_desired, q50 - q50_desired, q975 - q975_desired)

temps_stats = np.linspace(0.5, 4.5, 9) 


dist_params = {}
for impact in impacts.keys():
    dist_params[impact] = {}
        
    for t in temps_stats:
        
        q025_in = fit(t, *params[impact]['low'])
        q50_in = fit(t, *params[impact]['mid'])
        q975_in = fit(t, *params[impact]['high'])
    
        params_in = scipy.optimize.root(opt, [1, 1, 1], 
                    args=(q025_in, q50_in, q975_in)).x
            
        dist_params[impact][t] = params_in

#%%

percentiles = np.linspace(0.05, 0.95, 19)

percentiles = np.asarray([0.1, 0.5, 0.9])

linestyle_list = ['dotted', 'solid', 'dashed']

for impact in impacts.keys():
    
    for perc_i, percentile in enumerate(percentiles):
        vals = []
        for t in temps_stats:
            params_dist = dist_params[impact][t]
            
            vals.append(scipy.stats.skewnorm.ppf(percentile, 
                         params_dist[0], params_dist[1], params_dist[2]))
            
            
        params_percentile, _ = curve_fit(
            fit, temps_stats, vals)
             
    
        plt.plot(temps_plot, fit(temps_plot, *params_percentile), 
                 linestyle = linestyle_list[perc_i],
                 label=f'{impact} {100*percentile}', color=impacts[impact])

plt.xlabel(f'GMST cf pi (offset: {np.around(gmst_offset, decimals=3)}K)')
plt.ylabel('Mortality cf pi (solid), cf 2010 (dashed)')
plt.legend()
