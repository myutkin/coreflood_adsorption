#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 21:19:14 2021

@author: yutkinm
"""

from scipy.integrate import simps
import numpy as np
import matplotlib.pyplot as plt
import scipy as s

def calc_pv(tracer, d_vol, tracer_init, tracer_inj, cum_vol):
    cum_vol_bk = cum_vol
    cv = cum_vol_bk - d_vol
    tracer = (tracer - tracer_init)/(tracer_inj - tracer_init)
    tracer_init = (tracer_init - tracer_init)/(tracer_inj - tracer_init)
    tracer_inj = (tracer_inj - tracer_init)/(tracer_inj - tracer_init)
    A1 = list()
    A2 = list()
    for i in range(1,len(tracer)):
        I1 = simps(y=tracer[0:i], x=cv[0:i])
        A1.append(I1)
        print("A1 = {}".format(A1))
        I1 = simps(y = tracer[i:], x = cv[i:])
        I2 = tracer_inj*(cv[-1] - cv[i])
        A2.append(I2 - I1)
        print("A2 = {}".format(A2))
    A_diff = np.abs(np.array(A1) - np.array(A2))
    # print(A_diff)
    A_diff_min = np.min(np.abs(np.array(A1) - np.array(A2)))
    # this is a pv and dead volume together, now we need to get d_vol
    p_vol = cv[np.where(A_diff == A_diff_min)[0]]
    
    # This is for diagnostic, turn on when needed
    plt.figure()
    plt.scatter(cv,tracer)
    plt.vlines(p_vol, 0, 1)
    plt.hlines(1, 0, cv[-1])
    plt.hlines(0, 0, cv[-1])
#    plt.savefig('p_vol_diagnostic.png', dpi=72)
    plt.close('all')
    return p_vol

def make_opt(tracer_conc):
    def optimization(cum_vol, pe, d_vol, tracer_init, tracer_inj):
        tracer = tracer_conc
        global pp_vol
        p_vol = calc_pv(tracer, d_vol, tracer_init, tracer_inj, cum_vol) # calc_pv(tracer[-len(cum_vol):], cum_vol, d_vol)
        pp_vol = p_vol
        print(p_vol, d_vol, tracer_init, tracer_inj)
        x = (cum_vol - d_vol)/p_vol
        return (tracer_inj - tracer_init)*(0.5*(s.special.erfc((1-x)*s.sqrt(pe)/(2*s.sqrt(x)))) + 0.5*s.exp(pe)*s.special.erfc((1+x)*s.sqrt(pe)/(2*s.sqrt(x)))) + tracer_init
    return optimization

def plot_function(PV, pe, tracer_init, tracer_inj):
    x = PV
    return (tracer_inj - tracer_init)*(0.5*(s.special.erfc((1-x)*s.sqrt(pe)/(2*s.sqrt(x)))) + 0.5*s.exp(pe)*s.special.erfc((1+x)*s.sqrt(pe)/(2*s.sqrt(x)))) + tracer_init

