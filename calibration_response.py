from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TLorentzVector, TTree, TVector3
from ROOT import TH1D, TH2D, TFile, TTree, TColor, TCanvas, TLegend, TLatex, TLine, TMath, TGraphErrors, TF1, TProfile, TProfile2D
from ROOT import kBlack, kBlue, kRed, kYellow, kGreen, kGray, kOrange, kViolet
from ROOT import gStyle, gPad
from ROOT import gROOT
from ROOT import TStyle
from math import *
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from array import array
import os
import logging
import itertools
import fnmatch
import csv
import pandas as pd
import mplhep as hep

parser = ArgumentParser()
parser.add_argument('-i', '--inFiles', dest = 'inFiles', help="ntuples", nargs='+')
parser.add_argument('-p', '--particle', dest = 'particle', help="particle")
parser.add_argument("-o", "--outFile", dest = 'outFile',  help = "output file name")
args = parser.parse_args()

#load files
files = args.inFiles

#grab the particle type
particle = args.particle


# Energy binning, angular binning
EBins = array('d', (10.,50., 100., 150., 200., 250., 300., 350., 400., 450., 500., 550., 600., 650., 700., 750., 800.))#, 850., 900., 950., 1000.))
EBins_endcap = array('d', (10.,100., 200., 300., 400., 500., 600., 700., 800.))
EBins_lowE = array('d', (10.,15., 20., 25.,30.,35., 40.,45., 50.))
ThetaBins = np.linspace(0.225,2.91,30)
EBins_photon = array('d', (30.,50.,100.,150.,200.,250.,300.,400.,500.,650.,
    800.,1000., 1500., 2000., 2500., 3000., 3500., 4000., 4500., 5000.))
EBins_neutron = array('d', (25., 50., 75., 100., 150., 200., 250.))
phi_bins = np.linspace(-np.pi-7.5*np.pi/180., np.pi-7.5*np.pi/180.,25)
endcap_bins = np.linspace(0.225, 0.625, 5)
endcap_bins_2 = np.linspace(2.517, 2.91, 5)
tb_bins = np.linspace(0.625,2.517,20)
ebins_1 = np.array([0.175])
ebins_2 = np.array([2.96])



#####################################
# BINNED CALIBRATION                #
# double for loop, take the         #
# Etrue/Erec ratio for each theta/E #
# slice and apply to the 2D bins.   #
#####################################

corr_matrix = [] #for each theta bin, a list of avgs for each E bin
corr_matrix_tot = []
lowE_ratios=[]
lowE_thetas=[]
if particle == "neutron":
    EBins = EBins_neutron
if particle == "photon":
    EBins = EBins_photon
for tbin in range(0, len(ThetaBins)-1):
    ThMin = ThetaBins[tbin]
    ThMax = ThetaBins[tbin+1]
    Elist=[]
    for Ebin in range(0, len(EBins)-1):
        EMin = EBins[Ebin]
        EMax= EBins[Ebin+1]
        #if EMin > 30:
        #    continue
        ratios = []
        for filestr in files:
            file = TFile(filestr, "READ")
            tree = file.Get(f"{particle}_tree")
            for entry in tree:
                E_true = entry.E_truth
                theta = entry.theta_truth
                E_reco = entry.E
                theta_reco = entry.theta
                phi_reco = entry.phi
                if theta_reco < 0:
                    continue
                #use RECO values for Binning!!
                if EMin <= E_reco and E_reco < EMax:
                    if ThMin <= theta_reco and theta_reco < ThMax:
                        ratios.append(E_true / E_reco)
        if len(ratios) > 0:
            avg_ratio = np.average(ratios)
        else:
            #print("No entries in bin: ", tbin, " thmin: ", ThMin)
            avg_ratio = -1.
        Elist.append(avg_ratio)
    print("theta slice ", tbin, " done")
    print("Starting theta slice ",tbin+1,"...")
    print(Elist)
    corr_matrix.append(Elist)
#write calibration map to a csv file
outfile = args.outFile
with open(f'../ResponseMaps/{outfile}.csv', 'w', newline = '') as csvfile:
    mapwriter = csv.writer(csvfile)
    mapwriter.writerows(corr_matrix)
print(f"Created {outfile}.csv")

