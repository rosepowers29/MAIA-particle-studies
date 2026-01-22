from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TLorentzVector, TTree, TVector3
from ROOT import TH1F, TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TGraphErrors, TF1, TProfile, TProfile2D
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange, kViolet
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
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

parser = ArgumentParser()
parser.add_argument('-i', '--inFiles', dest = 'inFiles', help="ntuples", nargs='+')
parser.add_argument('-p', '--particle', dest = 'particle', help="particle")
parser.add_argument("-r", "--region", dest = 'region', default = 'e', help = "theta region")
parser.add_argument('-l', "--label", dest = 'label', default = 'NBv6', help = "outfile label (version, BIB or no BIB)")
args = parser.parse_args()

#load files
files = args.inFiles

#grab the region
region = args.region

#grab the particle type
particle = args.particle

#grab the label
label = args.label

# Energy binning (EPJC style), angular binning
EBins_photon = array('d', (30.,50.,100.,150.,200.,300.,400.,500.,650.,
    800.,1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500., 5000.))
EBins_neutron = array('d', (20.,40.,65.,100., 150., 200., 250.))

#ThetaBins = np.linspace(0.175,2.96,30)
NPhotonBins = np.linspace(0,50,50)


endcap_bins = np.linspace(0.225, 0.625, 5)
endcap_bins_2 = np.linspace(2.517, 2.91, 5)
tb_bins = np.linspace(0.625,2.517,20)
ThetaBins = np.concatenate((endcap_bins[:2], tb_bins, endcap_bins_2[1:]))

#declare arrays
e_arr = array('d')
th_arr = array('d')
th_err_arr=array('d')
sigma_arr = array('d')
sigma_arr_th = array('d')
e_err_arr = array('d')
sigma_err_arr = array('d')
sigma_err_arr_th = array('d')
mu_arr = array('d')
mu_err_arr = array('d')
mu_arr_th = array('d')
mu_arr_th_err = array('d')

################### RESOLUTION STUDY #################
cx = TCanvas("", "", 800, 600)
gStyle.SetOptStat(1)
Ecalfracs = []
EBins = None
if particle == "photon":
    EBins = EBins_photon
elif particle == "neutron":
    EBins = EBins_neutron
for Ebin in range(0, len(EBins)-1):
    #if EBins_v5[Ebin]>225:
    #    break
    if region == 'th':
        break
    EMin = EBins[Ebin]
    EMax = EBins[Ebin+1]
    proj_name = str(EBins[Ebin])+" GeV<E<"+str(EBins[Ebin+1])+" GeV"
    file_name = "Ereso"+str(Ebin)
    
    # adjust this binning as needed
    if EMin < 50.:
        lim = 1.
        bins =50
    elif 50<= EMin < 100:
        lim = 1.
        bins = 50
    elif 100<= EMin <= 400.:
        lim = .75
        bins = 50
    elif 400.< EMin < 1000:
        lim = .5
        bins = 50
    elif 1000. == EMin:
        lim = 0.75
        bins = 30
    else:
        lim = .5
        bins = 50

    h_my_proj_2 = TH1D(proj_name, proj_name, bins, -lim, lim)

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
            if E_reco < 10.:
                continue

            neg_lim = -10.
            pos_lim = 10.
            if (E_corr - E_truth)/E_truth > neg_lim and (E_corr - E_truth)/E_truth < pos_lim:
                h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
        file.Close()
    #now start fitting the gaussians

    lim = 1.
    if EMin<50000:
        gaussFit1 = TF1("gaussfit", "gaus", -lim, lim)
        gaussFit1.SetLineColor(kRed)
        gaussFit1.SetParameter(1, 0.)
        gaussFit1.SetParameter(2, h_my_proj_2.GetRMS())
        h_my_proj_2.Fit(gaussFit1, "E")
        gStyle.SetOptFit(0o111);
        x_axis = h_my_proj_2.GetXaxis()
        y_axis = h_my_proj_2.GetYaxis()
        x_axis.SetTitle("$\\Delta_{E}/E_{true}")
        y_axis.SetTitle("N Entries")
        h_my_proj_2.Draw("HIST")
        gaussFit1.Draw("LSAME")
        cx.Update()
        sigma = gaussFit1.GetParameter(2)
        sigma_err = gaussFit1.GetParError(2)
        bincenter = (EMax-EMin)/2
        e_arr.append(EMin+bincenter)
        e_err_arr.append(bincenter)
        sigma_arr.append(sigma)
        mu_arr.append(gaussFit1.GetParameter(1))
        mu_err_arr.append(gaussFit1.GetParError(1))
        print("SIGMA=",sigma)
        sigma_err_arr.append(sigma_err)
    

# theta scan
if region == 'th':

    for Tbin in range(0, len(ThetaBins)-1):
        ThMin = ThetaBins[Tbin]
        ThMax = ThetaBins[Tbin+1]
        proj_name = "E reso, "+str(np.round(ThetaBins[Tbin],2))+"<theta<"+str(np.round(ThetaBins[Tbin+1],2))
        file_name = "Ereso"+str(Tbin)
        lim = 0.5
        if ThMin > 0.577 and ThMin < 2.96:
            lim = 2.
        else:
            lim = 5.
        h_my_proj_2 = TH1D(proj_name, proj_name, 50,-lim,lim)

        for filestr in files:
            file = TFile(filestr, "READ")
            tree = file.Get(f"{particle}_tree")
            for entry in tree:
                E_truth = entry.E_truth
                theta = entry.theta_truth
                E_reco = entry.E 
                E_corr = entry.E
                theta_reco = entry.theta
                #get rid of negative theta values
                if theta_reco < 0:
                    continue
                # get rid of negative energy values, restrict range
                #only include the energy plateau
                #if E_truth < 600.:
                #    continue
                if E_reco < 10.:
                    continue
                neg_lim = -4.
                pos_lim = 4.
                if (E_corr - E_truth)/E_truth > neg_lim and (E_corr - E_truth)/E_truth < pos_lim:
                    h_my_proj_2.Fill((E_corr - E_truth)/E_truth)
            file.Close()
        #now start fitting the gaussians
        lim = 0.0
        if lim == 0.0:
            gaussFit1 = TF1("gaussfit", "gaus", -1., 1.)
            gaussFit1.SetLineColor(kRed)
            gaussFit1.SetParameter(1, 0.)
            gaussFit1.SetParameter(2, h_my_proj_2.GetRMS())
            h_my_proj_2.Fit(gaussFit1, "E")
            gStyle.SetOptFit(0o111);
            h_my_proj_2.Draw("HIST")
            gaussFit1.Draw("LSAME")
            cx.Update()
            #cx.SaveAs("slices_corr_full"+file_name+".root")
            sigma = gaussFit1.GetParameter(2)
            sigma_err = gaussFit1.GetParError(2)
            bincenter = ((ThMax-ThMin)/2)
            e_arr.append((ThMin + ThMax)/2)
            e_err_arr.append(bincenter)
            sigma_arr.append(sigma)
            mu_arr.append(gaussFit1.GetParameter(1))
            mu_err_arr.append(gaussFit1.GetParError(1))
            print("SIGMA=",sigma)
            sigma_err_arr.append(sigma_err)



# save arrays to a csv for easy access, plotting
# change filepath as needed
res_info = pd.DataFrame({'E': e_arr, 'sigma': sigma_arr, 'E_err': e_err_arr, 'sigma_err': sigma_err_arr})
res_info.to_csv(f"../ResArrays/{label}_resarrays_{particle}_{region}.csv")

