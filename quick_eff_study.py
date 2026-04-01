from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TH2D, TH1, TFile, TLorentzVector, TMath, TTree, TVector3, TEfficiency
from math import *
from argparse import ArgumentParser
from array import array
import os
import fnmatch
import uproot
import mplhep as hep
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument('-i', '--inFiles', dest = 'inFiles', help="ntuples", nargs='+')
parser.add_argument('-p', '--particle', dest = 'particle', help="particle")
parser.add_argument('-l', "--label", dest = 'label', default = 'NBv6', help = "outfile label (version, BIB or no BIB)")
args = parser.parse_args()

#load files
files = args.inFiles


#grab the particle type
particle = args.particle

#grab the label
label = args.label



arrBins_theta = np.linspace(0.175, TMath.Pi()-0.175, 30)#array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            #90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (15., 20., 25., 50., 100., 250., 500., 1000.))#, 2500., 5000.))
arrBins_E_lowrange = array('d', (15., 20., 25., 30., 35., 40., 45., 50., 75., 100., 125., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.))
arrBins_fakes = np.linspace(0,1000,250)
arrBins_dR = np.linspace(0,5.,100)
recobins = np.linspace(0, 10000, 100)

h_E_pass = TH1D('E_pass', 'E_pass', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_E_all = TH1D('E_all', 'E_all', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_theta_pass = TH1D('theta_pass', 'theta_pass', len(arrBins_theta)-1, arrBins_theta)
h_theta_all = TH1D('theta_all', 'theta_all', len(arrBins_theta)-1, arrBins_theta)


for filestr in files:
    file = TFile(filestr, "READ")
    tree = file.Get(f"{particle}_tree")
    for entry in tree:
        E_truth = entry.E_truth
        theta_truth = entry.theta_truth
    if theta_truth > 0:
        if theta_truth > 0.99 and theta_truth < 2.15:
            h_E_pass.Fill(E_truth)
        if E_truth > 10:
            h_theta_pass.Fill(theta_truth)
    if theta_truth > 0.175 and theta_truth < 2.96:
        h_E_all.Fill(E_truth)
    if E_truth > 20:
        h_theta_all.Fill(theta_truth)

outfile_name = f"quickeff_{args.particle}_{args.label}.root"
output_file = TFile(outfile_name, 'RECREATE')
h_theta_pass.Write()
h_theta_all.Write()
h_E_pass.Write()
h_E_all.Write()
output_file.Close()


