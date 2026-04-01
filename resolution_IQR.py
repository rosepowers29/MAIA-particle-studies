from ROOT import TFile, TTree
from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from math import *
import numpy as np
from argparse import ArgumentParser
from array import array
import os
import logging
import itertools
import fnmatch
import csv
import pandas as pd
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument('-i', '--inFiles', dest = 'inFiles', help="ntuples", nargs='+')
parser.add_argument('-p', '--particle', dest = 'particle', help="particle")
parser.add_argument("-r", "--region", dest = 'region', default = 'e', help = "theta region")
parser.add_argument('-l', "--label", dest = 'label', default = 'NBv6', help = "outfile label (version, BIB or no BIB)")
args = parser.parse_args()

#load files
files = args.inFiles
print(files)
#grab the region
region = args.region

#grab the particle type
particle = args.particle

#grab the label
label = args.label

# Energy binning
EBins_photon = array('d', (30.,50.,100.,200.,400.,600.,800.,
    1000., 1500., 2000., 2500., 3000., 4000., 5000.))
EBins_neutron = array('d', (20.,40.,65.,100., 150., 200., 250.))

#EBins_neutron = array('d', (10., 20., 30., 40., 50.))
#EBins_photon = array('d', (10., 20., 30., 40., 50.))
#ThetaBins = np.linspace(0.175,2.96,30)
NPhotonBins = np.linspace(0,50,50)


endcap_bins = np.linspace(0.225, 0.625, 5)
endcap_bins_2 = np.linspace(2.517, 2.91, 5)
tb_bins = np.linspace(0.625,2.517,20)
ThetaBins = np.concatenate((endcap_bins[:2], tb_bins, endcap_bins_2[1:]))

#declare arrays
e_arr = array('d')
e_err_arr = array('d')
iqr_res_arr = array('d')
iqr_res_err = array('d')
def prop_error(mu, sigma, mu_err, sigma_err):
    """
    Helper function to propagate the error on resolution measurement in quadrature
    """
    res = sigma / (mu+1) #scaled resolution
    res_err = res * np.sqrt((sigma_err / sigma)**2 + (mu_err / (mu+1))**2)
    return res_err

################### RESOLUTION STUDY #################
EBins = None
if particle == "photon":
    EBins = EBins_photon
elif particle == "neutron":
    EBins = EBins_neutron
bin_list = []
for Ebin in range(0, len(EBins)-1):
    if region == 'th':
        break
    EMin = EBins[Ebin]
    EMax = EBins[Ebin+1]
    print(EMin, "< E <", EMax)
    bin_list.clear()
    print("bin_list:", bin_list) #check to make sure it got cleared
    for filestr in files:
        file = TFile(filestr, "READ")
        tree = file.Get(f"{particle}_tree")
        for entry in tree:
            E_truth = entry.E_truth
            
            if E_truth < EMin or E_truth > EMax:
                continue
            theta = entry.theta_truth
            E_reco = entry.E
            E_corr = entry.E
            theta_reco = entry.theta
            phi_truth = entry.phi_truth
            phi_reco = entry.phi

            dR_reco_true = np.sqrt((phi_reco-phi_truth)**2+(theta_reco-theta)**2) 
            #restrict to endcap region (NO SOLENOID CONTACT)
    
            if region == 'e':
                if theta > 0.7 and theta <2.44 or theta < 0.275 or theta > 2.96:
                    continue
            #restrict to central-barrel region (minimal solenoid contact)
            elif region == 'b':
                if theta < .99 or theta > 2.15:
                    continue
            #restrict to transition region (maximal solenoid contact)
            elif region == 't':
                if theta < 0.7 or theta > 2.44 or 0.99 <theta< 2.15:
                    continue
            #restrict to barrel and trandition region
            elif region == 'tb':
                if theta < 0.6899 or theta > 2.45:
                    continue
            #all-inclusive
            elif region == 'a':
                if theta < 0.18 or theta > 2.96:
                   continue

            #get rid of any negative theta values
            if theta_reco < 0:
                continue
            # get rid of negative energy values, restrict range
            if E_reco < 20.:
                continue

            bin_list.append((E_corr - E_truth) / E_truth)
        file.Close()

    bincenter = (EMax-EMin)/2
    e_arr.append(EMin+bincenter)
    e_err_arr.append(bincenter)
        
    # calculated with 68% interquartile range
    print("len of bin_list after filling", len(bin_list))
    bin_arr = np.array(bin_list)
    med = np.median(bin_arr)
    iqr_68 = np.percentile(bin_arr, 84) - np.percentile(bin_arr, 16)
    
    iqr_res = abs(iqr_68 / (2 * (med + 1)))
    iqr_res_arr.append(iqr_res)

    # bootstrap error method (from Claude)
    n_boot = 10000
    boot_iqr68 = np.zeros(n_boot)
    boot_median = np.zeros(n_boot)
    boot_resolution = np.zeros(n_boot)

    for i in range(n_boot):
        sample = np.random.choice(bin_arr, n_boot, replace=True)
        q16, q84 = np.percentile(sample, [16, 84])
        med_ = np.median(sample)
        boot_iqr68[i] = (q84 - q16)
        boot_resolution[i] = boot_iqr68[i] /(2* (med + 1))  # already the full ratio

    err_resolution = boot_resolution.std()
    #err_iqr68     = boot_iqr68.std()
    #err_median    = boot_median.std()
    iqr_res_err.append(err_resolution)


    # plot individual hists
    plt.hist(bin_arr, bins = 50, color = "blue")
    plt.xlabel("$\Delta E / E_{{true}}$")
    plt.axvline(np.median(bin_arr), label = "Median", color = "red")
    plt.axvline(np.percentile(bin_arr, 16), label = "Percentile 16", color = "orange")
    plt.axvline(np.percentile(bin_arr, 84), label = "Percentile 84", color = "orange")
    plt.xlim(-0.6, 0.6)
    plt.legend()
    plt.savefig(f"distros/NB_{particle}_{region}_hist_{str(round(EMin, 0))}_{str(round(EMax, 0))}.png")
    plt.close()
res_info = pd.DataFrame({'E': e_arr, 'E_err': e_err_arr, 'iqr_res_err': iqr_res_err, 'iqr_res': iqr_res_arr})
res_info.to_csv(f"{label}_resoarrays_{particle}_{region}.csv")
print(f"Created file {label}_resoarrays_{particle}_{region}.csv")

