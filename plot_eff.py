from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D,TH1, TFile, TLorentzVector, TMath, TTree, TVector3, TEfficiency, TGraphAsymmErrors
from math import *
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import uproot
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

hasBIB=False

parser = ArgumentParser()
parser.add_argument('-n', '--NBhists', dest = 'NBhists', help="no-BIB ntuples", nargs='+')
parser.add_argument('-b', '--BIBhists', dest = 'BIBhists',help="BIB ntuples", nargs='+')
parser.add_argument('-v', '--var', dest = 'region', help="plotting variable")
parser.add_argument('-p', '--particle', dest = 'particle', help = "particle")
parser.add_argument('-l', '--label', dest = 'label', help = "label for output file")
args = parser.parse_args()

filedict = {"histofiles":[], "histofiles_B":[]}

if args.BIBhists is not None:
    for file in args.BIBhists:
        filedict["histofiles_B"].append(file)
    hasBIB=True
for file in args.NBhists:
    filedict["histofiles"].append(file)
print(args.region)
print(filedict["histofiles"])
print(filedict["histofiles_B"])
plotvar = args.region
particle = args.particle
label = args.label

def getErrorBars(eff_array, pass_slice, all_slice):
    """use the Clopper-Pearson method to obtain asymm error on the efficiency"""
    up_err_arr.clear()
    low_err_arr.clear()
    for i in range(len(eff_array)):
        N = pass_slice[i]
        k = all_slice[i]
        lower_error = eff_array[i] - TEfficiency.ClopperPearson(k, N,  0.683,  False)
        upper_error = TEfficiency.ClopperPearson(k, N,  0.683,  True) - eff_array[i] 
        up_err_arr.append(upper_error)
        low_err_arr.append(lower_error)
    asymm_errs = [np.asarray(low_err_arr[0:]), np.asarray(up_err_arr[0:])]
    return asymm_errs
filenum = 0
fig, ax = plt.subplots()
up_err_arr = []
low_err_arr = []
up_err_arr_B = []
low_err_arr_B = []
pass_slices = []
all_slices = []
E_slices = []
for filename in filedict["histofiles"]:
    if filenum >= len(filedict["histofiles"]):
        break
    file = uproot.open(filename)
    pass_histo = file[plotvar+'_pass;1']
    all_histo = file[plotvar+'_all;1']

    pass_histo_np = pass_histo.to_numpy()
    all_histo_np = all_histo.to_numpy()
    pass_histo_arr = np.asarray(pass_histo_np[0])
    all_histo_arr = np.asarray(all_histo_np[0])
    pass_slices.append(pass_histo_arr[0:])
    all_slices.append(all_histo_arr[0:])
    #all_histo_arr = all_histo_np[0]
    if filenum == len(filedict["histofiles"])-1:
        E_arr = pass_histo_np[1]
        E_slices.append(E_arr[0:])
    filenum = filenum+1
    file.close()

pass_slices_B = []
all_slices_B = []
E_slices_B = []
filenum = 0
if hasBIB:
    for filename in filedict["histofiles_B"]:
        if filenum >= len(filedict["histofiles_B"]):
            break
        file = uproot.open(filename)
        pass_histo = file[plotvar+'_pass;1']
        all_histo = file[plotvar+'_all;1']

        pass_histo_np = pass_histo.to_numpy()
        all_histo_np = all_histo.to_numpy()
        pass_histo_arr = np.asarray(pass_histo_np[0])
        all_histo_arr = np.asarray(all_histo_np[0])
        pass_slices_B.append(pass_histo_arr[0:])
        all_slices_B.append(all_histo_arr[0:])
        #all_histo_arr = all_histo_np[0]
        if filenum == len(filedict["histofiles_B"])-1:
            E_arr_B = pass_histo_np[1]
            E_slices_B.append(E_arr_B[0:])
        #print(E_slices_B)
        filenum = filenum+1
        print(filenum)
        file.close()
print(len(filedict["histofiles_B"]))


pass_slice = np.add.reduce(pass_slices)#pass_slice_0+pass_slice_1#+pass_slice_2
all_slice = np.add.reduce(all_slices)#all_slice_0+all_slice_1#+all_slice_2
efficiency_arr = pass_slice/all_slice

E_errs = []
E_arr = []
E_arr_B = [] 
E_errs_B = []

print(E_slices_B)
for i in range(len(E_slices[0])-1):
    E_arr.append((E_slices[0][i+1]+E_slices[0][i])/2)
    E_errs.append((E_slices[0][i+1]-E_slices[0][i])/2)
if hasBIB:
    for i in range(len(E_slices_B[0])-1):
            E_arr_B.append((E_slices_B[0][i+1]+E_slices_B[0][i])/2)
            E_errs_B.append((E_slices_B[0][i+1]-E_slices_B[0][i])/2)

#strip out the NaNs that arise from "empty" bins and the corresponding E values (i.e. for the 0-50 slice, all higher-E bins give us NaN eff)
E_arr = np.array(E_arr)
E_arr = E_arr[~np.isnan(efficiency_arr)]
clean_eff_arr = efficiency_arr[~np.isnan(efficiency_arr)]
E_errs = np.array(E_errs)
E_errs = E_errs[~np.isnan(efficiency_arr)]
if hasBIB:
    pass_slice_B = np.add.reduce(pass_slices_B)#np.concatenate(pass_slices_B, axis=0)#pass_slice_0_B + pass_slice_1_B + pass_slice_2_B
    all_slice_B = np.add.reduce(all_slices_B)#np.concatenate(all_slices_B, axis=0)#all_slice_0_B + all_slice_1_B + all_slice_2_B
    efficiency_arr_B = pass_slice_B / all_slice_B
    clean_eff_arr_B = efficiency_arr_B[~np.isnan(efficiency_arr_B)]
    E_arr_B = np.array(E_arr_B)
    E_arr_B = E_arr_B[~np.isnan(efficiency_arr_B)]
    E_errs_B = np.array(E_errs_B)
    E_errs_B = E_errs_B[~np.isnan(efficiency_arr_B)]
    asymm_errs_B = getErrorBars(clean_eff_arr_B, pass_slice_B, all_slice_B)
asymm_errs = getErrorBars(clean_eff_arr, pass_slice, all_slice) 

print(len(E_arr), ",", len(clean_eff_arr))
print(pass_slice)
print(all_slice)
plt.errorbar(E_arr, clean_eff_arr, xerr = E_errs[0:], yerr = asymm_errs,fmt=' ',
        color='blue', label = "Without BIB+IPP overlay")#label="BIB, jets, HF < 0.1, dR <= "+str(cutoffs[filenum]))
if hasBIB:
    print(len(E_arr[0:]), '\n', len(clean_eff_arr_B[0:]))
    plt.errorbar(E_arr_B, clean_eff_arr_B, xerr = E_errs_B, 
        yerr = asymm_errs_B, fmt = ' ', color = 'red', label="With BIB+IPP overlay") 
plt.axhline(1.,ls='--', lw=0.5)
ax.set_ylim(0.0,1.4)
maj_loc = 0
min_loc = 0
maj_form = 'null'
if plotvar == 'theta':
    ax.set_xlim(E_slices[0][0],E_slices[0][len(E_slices[0])-1])
    ax.set_xlabel(f"True {particle} $\\theta$ [rad]", loc='right')
    maj_loc = 0.5
    min_loc = 0.1
    maj_form = '{x:.1f}'
    plt.gcf().text(0.16, 0.7, f"20 GeV < True {particle} E < 250 GeV")
elif plotvar == 'E':
    #ax.set_xlim(E_slices[0][0],E_slices[0][len(E_slices[0])-1])
    ax.set_xlim(20,250)
    ax.set_xlabel(f"True {particle} Energy [GeV]", loc='right')
    maj_loc = 50
    min_loc = 10
    maj_form = '{x:.0f}'
    plt.gcf().text(0.16, 0.7, "0.175 < $\\theta$ < 2.96")
    #ax.set_xscale('log')
ax.set_ylabel(f"{particle} Matching Efficiency", loc='top')

hep.cms.label(exp = "Muon Collider", data = False, label = "with BIB+IPP (EU24 Lattice)", 
       rlabel='$MAIA$ Detector Concept', loc=2, italic=(1,0,0), pad=(0.0))
plt.gcf().text(0.16,0.74,"$\sqrt{s}=10$ TeV")

ax.legend(frameon=False)

ax.tick_params(bottom=True, top=True, left=True, right=True, which='both', direction='in')
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)
ax.tick_params(axis= 'y',which='major',length = 10)
ax.tick_params(axis='y',which='minor', length=5)

#ax.xaxis.set_major_locator(MultipleLocator(maj_loc))
#ax.xaxis.set_major_formatter(maj_form)
#ax.xaxis.set_minor_locator(MultipleLocator(min_loc))

ax.yaxis.set_major_locator(MultipleLocator(0.1))
ax.yaxis.set_major_formatter('{x:.1f}')
ax.yaxis.set_minor_locator(MultipleLocator(0.02))
plt.savefig(f"plots/{particle}Eff_{label}_maia_"+plotvar+".pdf")
plt.close()
print(f"created file plots/{particle}Eff_{label}_maia_"+plotvar+".pdf")
