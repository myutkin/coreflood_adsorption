#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 17:07:25 2019

@author: yutkinm
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 14:41:40 2019

@author: yutkinm
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.stats.distributions import  t
import pandas as pd
import functions
import matplotlib

matplotlib.font_manager._rebuild()

import os
os.environ["PATH"] += os.pathsep + '/Users/maxim/bin'

print(os.environ['PATH']) 

#%%

#%% auxilary finctions and definitions

plt.rcParams['font.size'] = 30
plt.rcParams['figure.figsize'] = (16,7);
plt.rcParams['figure.figsize'] = (15, 6.5625);
plt.rcParams['figure.figsize'] = (20, 10);
plt.rcParams['lines.linewidth'] = 3;
plt.rcParams['lines.markersize'] = 7;
plt.rcParams['axes.labelsize'] = 30;
plt.rcParams['xtick.labelsize'] = 30;
plt.rcParams['ytick.labelsize'] = 30;
plt.rcParams['legend.fontsize']=17;
plt.rcParams["legend.handlelength"] = 3;
plt.rcParams["legend.borderpad"] = 0.3;
plt.rcParams["legend.labelspacing"] = 0.4;

plt.rcParams['font.family'] = 'serif' 
plt.rcParams['font.serif'] = 'Computer Modern'


plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.preamble'] = r''' \usepackage[version=4]{mhchem} 
 \usepackage{siunitx} 
 \usepackage{relsize}'''

#%%

# Load your data on ceoncentration and incremental volume

filename = 'SDS_adsorption_IL-HP-4_Na Ca S Ksenia 25.8.2021'
DATA = pd.read_excel(filename + '.xlsx')

# d_vol = 9 # mL, measured experimentally on the coreflood
#

CV_list = len(DATA['Ca_mmol/L'].dropna())*[2]

cum_vol1 = np.cumsum(CV_list) # S data


CV_list.insert(0,1)
CV_list.pop(-1)

cum_vol2 = np.cumsum(CV_list) # Cl data


A = 11.28 # cm2
L = 6.94 # cm

V_tot = A * L # total volume of the plug!!, ml
#
#volume_step = np.array(DATA['Volume, ml'].values.tolist())
#volume_step = volume_step.astype(float)    
#
#cum_vol = np.cumsum(volume_step)
cum_vol = cum_vol2
cum_vol_glob = cum_vol

# load concentration for each volume step, concentration can be in any units 
# (except dimensionless), make sure to provide initial guess for optimization 
# algorithm in the same units as here (see below)

tracer = np.array(DATA['Cl, ppm'].dropna().values.tolist()) 
tracer = tracer.astype(int)   

S_exp = np.array(DATA['S_mmol/L'].dropna().values.tolist()) 

# p0 is initial guess for Peclet, dead volume, initial tracer concentration, 
# and injected tracer concentration
    
# bounds are [lower] and [upper] bounds for the parameters in the same order

result = curve_fit(functions.make_opt(tracer[5:]), cum_vol[5:], tracer[5:], p0=(5, 9, 3700, 380), 
                   bounds=([2, 7, 3300, 300], [20, 10, 4000, 500]))

pe = result[0][0]
d_vol = result[0][1]
tracer_init = result[0][2]
tracer_inj = result[0][3]

perr = np.sqrt(np.diag(result[1]))

alpha = 0.05 # 95% confidence interval = 100*(1-alpha)

n = len(tracer[5:])    # number of data points
p = 4 # number of parameters
dof = max(0, n - p) # number of degrees of freedom
# student-t value for the dof and confidence level
tval = t.ppf(1.0-alpha/2., dof) 

p_vol = functions.pp_vol
PV1 = (cum_vol1-d_vol)/p_vol
PV2 = (cum_vol2-d_vol)/p_vol


# this is rescaling since the model expects an input from 0 to 1

tracer_norm = (tracer - tracer_init)/(tracer_inj - tracer_init)

#%% plot
plt.close('all')

fig = plt.figure()

ax1 = fig.add_subplot(1,2,1)

ax1.yaxis.tick_left()
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(axis='y', which='both', direction='in') 


points1 = ax1.scatter(PV2, tracer_norm, ec='blue', fc='none', s=150)

line1 = ax1.plot(np.arange(0.01, 5, 0.01), 
         ((functions.plot_function(np.arange(0.01, 5, 0.01), pe, 
                         tracer_init, tracer_inj)) - tracer_init)/(tracer_inj - tracer_init),
         c='blue')
         
points2 = ax1.scatter(PV1, S_exp/16.1875, ec='red', fc='none', s=150)    
         
ax1.grid(linewidth=0.5, linestyle = '--')

textstr = '\n'.join((
        'IL-HP-4',
        r'$Pe = %.1f \pm %.1f$' % (pe, perr[0], ),
#        r'$V_d = %.1f \pm %.1f ~ \mathrm{mL}$' % (d_vol, perr[1]*tval, ),
        r'$\varphi_F = %.2f \pm %.3f$' % (functions.pp_vol/V_tot, 2/V_tot),
#        r'$Tracer_{inj} = %.1f \pm %.1f $' % (tracer_inj, perr[3]*tval, ),
#        r'$Tracer_{init} = %.1f \pm %.1f $' % (tracer_init, perr[2]*tval, )
        ))
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
#ax1.text(0.5, 0.2, textstr, fontsize=24, bbox=props,transform=fig.transFigure)
ax1.text(2, 0.3, textstr, fontsize=28, bbox=props)

ax1.set_xlabel(r'Dimensionless time, $ut/\varphi L$')
ax1.set_ylabel(r'Dimensionless concentration, $\widetilde{C}$')
ax1.set_xlim([-0.1, 5.1])
ax1.set_xticks([0, 1, 2, 3, 4, 5])

ax1.set_ylim([-0.01, 1.01])
ax1.set_yticks([0, 0.1, 0.3, 0.5, 0.7, 0.9, 1])

plt.show('fig')
fig.tight_layout()
fig.savefig(filename + '_tracer.png', dpi=300)

# http://kitchingroup.cheme.cmu.edu/blog/2013/02/12/Nonlinear-curve-fitting-with-parameter-confidence-intervals/