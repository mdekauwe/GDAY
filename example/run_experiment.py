#!/usr/bin/env python

""" run GDAY, plot LAI, NPP, transpiration """

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import datetime as dt

__author__  = "Martin De Kauwe"
__version__ = "1.0 (11.02.2014)"
__email__   = "mdekauwe@gmail.com"

def main(site, treatment):
    
    # run new simulations
    os.system("example.py")

    # load data
    amb = read_data("outputs/D1GDAY%s%s.csv" % (site, treatment.upper()))
    
    plt.rcParams['figure.subplot.hspace'] = 0.15
    plt.rcParams['figure.subplot.wspace'] = 0.15
    plt.rcParams['font.size'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['xtick.labelsize'] = 10.0
    plt.rcParams['ytick.labelsize'] = 10.0
    plt.rcParams['axes.labelsize'] = 10.0
    plt.rcParams['font.family'] = "sans-serif"
    plt.rcParams['font.style'] = "normal"
    plt.rcParams['font.serif'] = "Helvetica"
    fig = plt.figure()

    #years = np.unique(amb.YEAR).values
    years = np.unique(amb.YEAR)
    
    ax1 = fig.add_subplot(311)
   
    ax1.plot(years, amb.groupby("YEAR").NPP.sum(), "b-", label="Amb")
    
    obs_npp = np.array([1010.363435,1052.484494,930.9501267,1132.113908,\
                        1204.210228,1020.08473,662.2080009,953.6588582,\
                        1038.468081,909.183305,1220.6568,965.8984])
    ax1.plot(years, obs_npp, "ro", label="Obs")
    ax1.legend(numpoints=1, loc="best")
    ax1.set_ylabel("NPP (g m$^{2}$ d$^{-1}$)")
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_ylim(0., 1500)
    ax1.set_xlim(1995, 2008)
    
    ax2 = fig.add_subplot(312)
    ax2.plot(years, amb.groupby("YEAR").LAI.max(), "b-", label="Amb")
    #ax2.legend(numpoints=1, loc="best")
    ax2.set_ylabel("LAI (m$^{2}$ m$^{-2}$)")
    ax2.set_xlabel("Year")
    ax2.set_ylim(0., 5)
    
    ax3 = fig.add_subplot(313)
    ax3.plot(years, amb.groupby("YEAR").T.sum(), "b-", label="Amb")
    ax3.legend(numpoints=1, loc="best")
    ax3.set_ylabel("Transpiration (mm d$^{-1}$)")
    ax3.set_xlabel("Year")
    ax3.set_ylim(0., 500)
    plt.show()
    #fig.savefig("example_plot.png", dpi=100)
    
    
def date_converter(*args): 
    return dt.datetime.strptime(str(int(float(args[0]))) + " " +\
                                str(int(float(args[1]))), '%Y %j')

def read_data(fname):
    df = pd.read_csv(fname, parse_dates=[[0,1]], index_col=0, sep=",", 
                     keep_date_col=True, date_parser=date_converter, 
                     na_values=["-999.9"], skiprows=3)
    return df
    
    


if __name__ == "__main__":
    
    site = "DUKE"
    treatment="amb"
    main(site, treatment)
