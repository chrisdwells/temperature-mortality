import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import copy
import scipy.stats
from matplotlib.lines import Line2D
import pickle

figdir = '../figures'
datadir = '../data'

def fit(x, a, beta_t, beta_t2):
    yhat = a + beta_t*x + beta_t2*x**2
    return yhat

impacts = {
    'hot':'red', 
    'cold':'blue', 
    # 'net':'black',
    }

percs = ['low', 'mid', 'high']

gmst = f'{datadir}/AR6_GMST.csv'
gmst_in = pd.read_csv(gmst)

gmst_2000_2019 = np.mean(gmst_in[(2000 <= gmst_in['year']) & 
                     (gmst_in['year']<= 2019)]['gmst'])

gmst_offset = gmst_2000_2019
# gmst_1970_1999 = np.mean(gmst_in[(1970 <= gmst_in['year']) & 
#                      (gmst_in['year']<= 1999)]['gmst'])

# gmst_offset = gmst_2000_2019 - gmst_1970_1999

#%%

temps_plot = np.linspace(0, 5, 100)

params = {}

for impact in impacts.keys():
    params[impact] = {}
    
    for perc in percs:
            
        mort_file = f'{datadir}/Bressler_fig3_{impact}_{perc}.csv'
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

plt.tight_layout()
plt.savefig(f'{figdir}/temp_mort.png', dpi=100)
# plt.clf()    

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

params_percentiles = {}

percentiles = np.linspace(0.025, 0.975, 11)

for impact in impacts.keys():
    params_percentiles[impact] = {}
    
    for perc_i, percentile in enumerate(percentiles):
        vals = []
        for t in temps_stats:
            params_dist = dist_params[impact][t]
            
            vals.append(scipy.stats.skewnorm.ppf(percentile, 
                         params_dist[0], params_dist[1], params_dist[2]))
            
            
        params_percentile, _ = curve_fit(
            fit, temps_stats, vals)
        
        params_percentiles[impact][percentile] = params_percentile
        
#%%


for impact in impacts.keys():
    
    for perc in percs:
            
        mort_file = f'{datadir}/Bressler_fig3_{impact}_{perc}.csv'
        df_in = pd.read_csv(mort_file)
        
        temps = df_in['Temperature'].values
        morts = df_in['Mortality'].values
    
        temps_offset = temps + gmst_offset
    
        plt.scatter(temps_offset, morts, color=impacts[impact])
    
    
        params_in, _ = curve_fit(
            fit, temps_offset, morts)
            
        plt.plot(temps_plot, fit(temps_plot, *params_in), color=impacts[impact],
                 linestyle='--')
        

    plt.xlabel(f'GMST cf pi (offset: {np.around(gmst_offset, decimals=3)}K)')
    plt.ylabel('Mortality cf pi (solid), cf 2001-20 (dashed)')
    # plt.legend()

    
    for perc_i, percentile in enumerate(percentiles):

        params_in = params_percentiles[impact][percentile]             
    
        plt.plot(temps_plot, fit(temps_plot, *params_in), 
                 # linestyle = linestyle_list[perc_i],
                 label=f'{impact} {100*percentile}', color=impacts[impact])

    # plt.xlabel(f'GMST cf pi (offset: {np.around(gmst_offset, decimals=3)}K)')
    # plt.ylabel('Mortality cf pi (solid), cf 2010 (dashed)')
    # plt.legend()


handles = []
for impact in impacts.keys():
    handles.append(Line2D([0], [0], label=impact, marker='.', markersize=10, 
         markeredgecolor=impacts[impact], markerfacecolor=impacts[impact], linestyle=''))


handles.append(Line2D([0], [0], label='Rebased, 2.5-97.5', color='grey'))
handles.append(Line2D([0], [0], linestyle='--', label='Quadratic, 2.5-97.5', color='grey'))

plt.legend(handles=handles)

    
plt.tight_layout()
plt.savefig(f'{figdir}/temp_mort_percentiles.png', dpi=100)
# plt.clf()     


#%%

with open(f'{datadir}/outputs/params.pickle', 'wb') as handle:
    pickle.dump(params_percentiles, handle, protocol=pickle.HIGHEST_PROTOCOL)

#%%

output_dict = {}
output_dict['Percentiles'] = percentiles
output_dict['Demographics.Cold mortality sensitivity to T[1]'
            ] = np.full(percentiles.shape, np.nan)
output_dict['Demographics.Cold mortality sensitivity to T2[1]'
            ] = np.full(percentiles.shape, np.nan)
output_dict['Demographics.Hot mortality sensitivity to T[1]'
            ] = np.full(percentiles.shape, np.nan)
output_dict['Demographics.Hot mortality sensitivity to T2[1]'
            ] = np.full(percentiles.shape, np.nan)

for perc_i, percentile in enumerate(percentiles):
    output_dict['Demographics.Cold mortality sensitivity to T[1]'
                ][perc_i] = params_percentiles['cold'][percentile][1]*0.01 # from % to fraction
    
    output_dict['Demographics.Cold mortality sensitivity to T2[1]'
                ][perc_i] = params_percentiles['cold'][percentile][2]*0.01
    
    output_dict['Demographics.Hot mortality sensitivity to T[1]'
                ][perc_i] = params_percentiles['hot'][percentile][1]*0.01
    
    output_dict['Demographics.Hot mortality sensitivity to T2[1]'
                ][perc_i] = params_percentiles['hot'][percentile][2]*0.01


with open(f'{datadir}/outputs/mortality_output_dict.pickle', 'wb') as handle:
    pickle.dump(output_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


