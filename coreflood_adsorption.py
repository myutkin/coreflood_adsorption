#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:09:31 2021

@author: maxim
"""

#%%

import numpy as np
from subprocess import Popen

# import os
# cwd = os.getcwd()
# os.environ["PATH"] += os.pathsep + '/Users/maxim/bin'
# print(os.environ['PATH']) 

import pandas as pd
import matplotlib.pyplot as plt

phreeqc_bin = '/home/z/bin/phreeqc'
phreeqc_database = '/home/z/share/doc/phreeqc/database/phreeqc.dat'

#%%  plotting settings

plt.rcParams['font.size'] = 30
plt.rcParams['figure.figsize'] = (16,7);
plt.rcParams['figure.figsize'] = (15, 6.5625);
plt.rcParams['figure.figsize'] = (10, 10);
plt.rcParams['lines.linewidth'] = 3;
plt.rcParams['lines.markersize'] = 7;
plt.rcParams['axes.labelsize'] = 30;
plt.rcParams['xtick.labelsize'] = 30;
plt.rcParams['ytick.labelsize'] = 30;
plt.rcParams['legend.fontsize']=26;
plt.rcParams["legend.handlelength"] = 3;
plt.rcParams["legend.borderpad"] = 0.3;
plt.rcParams["legend.labelspacing"] = 0.4;

plt.rcParams['font.family'] = 'serif' 
plt.rcParams['font.serif'] = 'Computer Modern'


plt.rcParams['text.usetex']=True
plt.rcParams['text.latex.preamble'] = r''' \usepackage[version=4]{mhchem} 
 \usepackage{siunitx} 
 \usepackage{relsize}'''


#%% exp params
# sample name IL-HP-4

Q = 8.33333e-9 # m3/s, this is 0.50 mL/min
phi = 0.147 # by tracer test
L = 69.41e-3 # m, avg of 5
D = 37.9e-3 # m, avg of 5    
A = np.pi*D**2/4
u = Q/A # m/s

Kl = 32e-3 # Da
a_v = np.sqrt(36*phi**3/(Kl*9.869e-13*150*(1-phi)**2)) # well-known Karmen-Cozeney
D_g = 6/a_v # m, approximation for spheres

N_av = 6.02e23 # 1/mol

PV = 11.3e-6 # A*L*phi # m3
# This is an ion exchange experiment
# but we first assume chloride does not participate and 
# only SDS and carbonate exchange,
# we still treat Cl as a tracer

Na_init = 0.1 # M
Na_inj= 0.024 # M
Cl_init = Na_init
Cl_inj = 0.01

# phreeqc_input = phreeqc_input_common + '\n' + phreeqc_input_loading + '\n' + phreeqc_input_export

#%% phreeqc input

PC_input = '''
# calculation of the initial and injected solutions


SOLUTION_MASTER_SPECIES
Sds   Sds-   1  Sds  265

SOLUTION_SPECIES

Sds- = Sds-
log_k 0

SOLUTION 10000 initial solution in the core
units mol/L
pH 7 charge
Na {Na_init}
Cl {Cl_init}
# K 1e-3
# C(4) 0.5e-3


# first equilibrate with calcite-air
EQUILIBRIUM_PHASES 10000
Calcite 0 10
CO2(g) -3.35 10

# then take above and equilibrate with calcite only
# this is our composition inside!
EQUILIBRIUM_PHASES 10001
Calcite 0 10

SAVE SOLUTION 10000
COPY SOLUTION 10000 1-{N}

SOLUTION 0 Injected solution
pH        7 charge
units     mol/L

Na {Na_inj}
Cl {Cl_inj}
Sds 0.014

EQUILIBRIUM_PHASES 0
Calcite 0 10
CO2(g) -3.5 10

#SAVE SOLUTION 20000
#COPY SOLUTION 20000 0

#print
#-reset false
END
              EXCHANGE_MASTER_SPECIES
              Z Z-
              Y Y+
              EXCHANGE_SPECIES
              Z- = Z- ; log_k 0
              Y+ = Y+ ; log_k 0

              Z- + Na+ = NaZ
              log_k 0

              2Z- + Ca+2 = CaZ2
              log_k {log_k}
              
              Y+ + HCO3- = YHCO3
              log_k 0
              
              #2Y+ + CO3-2 = Y2CO3
              #log_k 0
              
              Y+ + Sds- = YSds
              log_k {log_k1}
              

              EXCHANGE 1-{N}
              -equilibrate 1
               Z   {CEC}
               Y   {CEC1}

           EQUILIBRIUM_PHASES 1-{N}
       		Calcite 0 2
         
           TRANSPORT

            -boundary_conditions	  flux   flux #
            -dispersivities		    {N}*{disp}
            -correct_disp		      false
            -diffusion_coefficient	 2e-9
            -cells {N}
            -shifts {shifts}
            -length {N}*{LN}
            -punch_cells {N}
            -punch_frequency 	   1
            -print_cells 1

            PRINT
            # -reset false
            -user_print            true
            -selected_output       true
            #-status                false
            -warnings              100


            USER_PUNCH 1
            -headings mCa mNa mCl mHCO3 mCO3 nH mOH
            -start
            1 m_h2o = 18e-3
            10 b0 = 1/m_h2o
            20 c0 = RHO * b0/(1 + MOL("Na+")*23e-3 + MOL("Ca+2")*40e-3 + MOL("Cl-")*35.5e-3 + MOL("CO3-2")*60e-3 + MOL("HCO3-")*61e-3 + MOL("OH-")*17e-3 + MOL("H+")*1e-3)
            30 mCa = c0 * m_h2o * MOL("Ca+2") * ACT("H2O")
            40 PUNCH mCa
            50 mNa = c0 * m_h2o * MOL("Na+") * ACT("H2O")
            60 PUNCH mNa
            70 mCl = c0 * m_h2o * MOL("Cl-") * ACT("H2O")
            80 PUNCH mCl
            90 mHCO3 =  c0 * m_h2o * MOL("HCO3-") * ACT("H2O")
            100 PUNCH mHCO3
            110 mCO3 =  c0 * m_h2o * MOL("CO3-2") * ACT("H2O")
            120 PUNCH mCO3
            130 mH = c0 * m_h2o * MOL("H+") * ACT("H2O")
            140 PUNCH mH
            150 mOH = c0 * m_h2o * MOL("OH-") * ACT("H2O")
            160 PUNCH mOH
            -end


            SELECTED_OUTPUT 1
            -file {filename}
            -high_precision       true
            -reset                false
            -distance             true
            -step                 true
            -pH                   true
            -totals               Na  Ca Cl Sds
            END            
            
''' 

#%% compose input vars, create pc input file, and run pc on it

N = 50 # number of cells
LN = L/N
#peclet = 100 # for testing to speed up calcs
peclet = 3 # see tracer calcs and graph

disp = L/peclet
CEC = 0.03 # 5*a_v*1e18/N_av/phi/1000 # mol/L, cation CEC
CEC1 = 0.05 # mol/L from tracer data #0.3*a_v*1e18/N_av/phi/1000 # mol/L, anion CEC
log_k = 4 # cation exch const
log_k1 = -1.25 # anion exchange const
shifts = 4*N # 4 pv
filename = 'IL-HP-4'
filename_ex = './' + filename + '.pco'

input_vars = {
       "N":N,
       "disp":disp,
       "LN":LN,
       "log_k":log_k,
       "CEC":CEC,
       "filename":filename_ex,
       "Na_inj":Na_inj,
       "Cl_inj":Cl_inj,
       "Na_init":Na_init,
       "Cl_init":Cl_init,
       "shifts":shifts,
       "log_k1":log_k1,
       "CEC1":CEC1
       }

with open(filename + '.pci', "w") as text_file:
   text_file.write(PC_input.format(**input_vars))


p = Popen([phreeqc_bin, './' + filename + '.pci', "OUTPUT.tmp", phreeqc_database])
p.wait()

#%% import data: exp and calc

exp_data = pd.read_excel('SDS_adsorption_IL-HP-4_Na Ca S Ksenia 25.8.2021.xlsx', header=0, skiprows=list(range(1,24)))

d_vol = 7 # mL, measured experimentally but also fitted to tracer exp, which is preferred

CV_list = len(exp_data['Ca_mmol/L'].dropna())*[2]

cum_vol1 = np.cumsum(CV_list) - d_vol # for ICP data: Ca, Na, Cl
CV_list.insert(0,1)
CV_list.pop(-1)

cum_vol2 = np.cumsum(CV_list) - d_vol # for Cl data

pc_data = pd.read_csv(filename_ex, sep="\t", 
                                    skipinitialspace=True, header=0,
                                    skiprows = lambda x: x in [1, 2] )
pc_data.columns = pc_data.columns.str.strip()

#%% plot the data

PV_calc = (pc_data['step'] + 0.5)/N
Na_calc = pc_data['mNa']
Ca_calc = pc_data['mCa']
Cl_calc = pc_data['mCl']
S_calc = pc_data['Sds']

PV_exp1 = cum_vol1/(PV * 1e6)
PV_exp2 = cum_vol2/(PV * 1e6)

Na_exp = exp_data['Na_mol/L'].dropna() # mol/L
Ca_exp = exp_data['Ca_mmol/L'].dropna() # mmol/L
Cl_exp = exp_data['Cl, ppm'].dropna()/35.5/1000 # mol/L
S_exp = exp_data['S_mmol/L'].dropna() # mmol/L

plt.close('all')
f1, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.tick_params(axis='y', which='both', direction='in', pad=int(7)) 
ax1.tick_params(axis='x', which='both', direction='in', top='on') 
ax2.tick_params(axis='y', which='both', direction='in', pad=int(7)) 

l1 = ax1.plot(PV_calc, Na_calc, c='blue', linewidth=2, linestyle='-.', label=r'\ce{Na^+}$_{\text{model}}$')
l2 = ax1.scatter(PV_exp1, Na_exp, marker='s', ec='blue', fc='none', linewidth=3, s=150, label=r'\ce{Na^+}$_{\text{exp}}$')

l5 = ax1.plot(PV_calc, Cl_calc, c='green', linewidth=2, linestyle='--', label=r'\ce{Cl^-}$_{\text{model}}$')
l6 = ax1.scatter(PV_exp2, Cl_exp, marker='>', ec='green', fc='none', linewidth=3, s=150, label=r'\ce{Cl^-}$_{\text{exp}}$')


ax1.set_ylabel(r'\ce{Cl^-}, \ce{Na^+} concentration, M')
ax1.set_xlabel(r'Dimensionless time, $ut/\varphi L$')

#ax1.set_ylim([0, 0.12])
#ax1.set_xlim([-0.01, 4])
ax1.grid(linewidth=0.5, linestyle=':')

l3 = ax2.plot(PV_calc, Ca_calc*1000 * 10, c='red', linewidth=2, linestyle=':', label=r'\ce{Ca^2+}$_{\text{model}}$')
l4 = ax2.scatter(PV_exp1, Ca_exp * 10, marker='o', ec='red', fc='none', linewidth=3, s=150, label=r'\ce{Ca^2+}$_{\text{exp}}$')

l7 = ax2.plot(PV_calc, S_calc*1000, c='black', linewidth=2, label=r'SDS$_{\text{model}}$')
l8 = ax2.scatter(PV_exp1, S_exp, marker='^', ec='black', fc='none', linewidth=3, s=150, label=r'SDS$_{\text{exp}}$')

ax2.set_ylabel(r'\ce{Ca^2+}, \ce{R-SO4^-} concentration, mM')
#ax2.set_ylim([6e-1, 1.2])

lns = [l1[0], l2, l3[0], l4, l5[0], l6, l7[0], l8]
labs = [i.get_label() for i in lns]
ax1.legend(lns, labs, loc='lower left', bbox_to_anchor=(0.05,1,0,0), 
           frameon=True, labelspacing=0.25, borderpad=0.25, numpoints=1, ncol=2)

#f1.tight_layout()
f1.savefig('coreflood.png', dpi=300)
plt.show()