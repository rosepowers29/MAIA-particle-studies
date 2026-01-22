from pyLCIO import IOIMPL
from pyLCIO import EVENT, UTIL
from ROOT import TH1D, TH2D, TH1, TFile, TLorentzVector, TMath, TTree, TVector3, TEfficiency
from math import *
from optparse import OptionParser
from array import array
import os
import fnmatch
import uproot
import mplhep as hep
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#-------------------------------------
def find_dr(phi1, theta1, phi2, theta2):
    dphi = abs(phi1-phi2)
    dtheta = abs(theta1-theta2)
    dr = np.sqrt(dphi**2+dtheta**2)
    return(dr)
#########################
parser = OptionParser()
parser.add_option('-i', '--inFile', help='--inFile Output_REC.slcio',
                  type=str, default='Output_REC.slcio')
parser.add_option('-o', '--outFile', help='--outFile ntup_photons.root',
                  type=str, default='ntup_photons.root')
(options, args) = parser.parse_args()

arrBins_theta = np.linspace(0.175, TMath.Pi()-0.175, 30)#array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            #90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (15., 20., 25., 50., 100., 250., 500., 1000.))#, 2500., 5000.))
arrBins_E_lowrange = array('d', (15., 20., 25., 30., 35., 40., 45., 50., 75., 100., 125., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.))
arrBins_fakes = np.linspace(0,1000,250)
arrBins_dR = np.linspace(0,5.,100)

h_E_pass = TH1D('E_pass', 'E_pass', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_E_all = TH1D('E_all', 'E_all', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_theta_pass = TH1D('theta_pass', 'theta_pass', len(arrBins_theta)-1, arrBins_theta)
h_theta_all = TH1D('theta_all', 'theta_all', len(arrBins_theta)-1, arrBins_theta)
# declare histograms
h_truth_E = TH1D('truth_E', 'truth_E', len(arrBins_E)-1, arrBins_E)
h_truth_theta = TH1D('truth_theta', 'truth_theta', len(arrBins_theta)-1, arrBins_theta)

h_truth_pT = TH1D('truth_pT', 'truth_pT', len(arrBins_E)-1, arrBins_E) #GeV: true generator pT

h_reco_neutron_E = TH1D('reco_photon_E', 'reco_photon_E', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)

h_hit_dR = TH1D('dR', 'dR', len(arrBins_dR)-1, arrBins_dR)
h_min_hit_dR = TH1D('mindR', 'mindR', len(arrBins_dR)-1, arrBins_dR)
h_hit_r = TH1D('hit_r', 'hit_r', len(arrBins_dR)-1, arrBins_dR)
histos_list = [h_truth_E, h_truth_theta, h_truth_pT, h_reco_neutron_E, h_hit_dR,
        h_min_hit_dR, h_hit_r,
        h_E_pass, h_E_all,
        h_theta_pass, h_theta_all]
# Histo list for writing to outputs

for histo in histos_list:
    histo.SetDirectory(0)
n_pass=0
n_all=0
####################################
neutron_tree = TTree("neutron_tree", "neutron_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
ecal_frac = array('d', [0])
neutron_tree.Branch("E",  E,  'var/D')
neutron_tree.Branch("phi", phi, 'var/D')
neutron_tree.Branch("theta", theta, 'var/D')
neutron_tree.Branch("E_truth",  E_truth,  'var/D')
neutron_tree.Branch("theta_truth", theta_truth, 'var/D')
neutron_tree.Branch("phi_truth", phi_truth, 'var/D')
neutron_tree.Branch("ecal_frac", ecal_frac, 'var/D')


to_process = []

if os.path.isdir(options.inFile):
    for r, d, f in os.walk(options.inFile):
        for file in f:
            to_process.append(os.path.join(r, file))
else:
    to_process.append(options.inFile)

filenum=0
nevts_proc = 0
for file in to_process:
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    try:
        reader.open(file)
    except Exception:
        #let it skip the bad files without breaking
        print("skipping file", filenum)
        continue
    filenum=filenum+1
    # loop over all events in the file
    for ievt, event in enumerate(reader):
        #if filenum % 10 == 0:
        print(" ")
        print("File "+str(filenum))
        print("Processing event " + str(ievt))



        #get the truth particle, always the first in the collection
        mcpCollection = event.getCollection('MCParticle')
        hasBadNuclei = False
        for mcp in mcpCollection:
            if mcp.getPDG() > 100000: #nuclei have 10-digit pdgid
                if mcp.getEnergy() > 4: #greater than 4 GeV
                    hasBadNuclei = True
                    break

        if hasBadNuclei == True:
            continue
        trueNeutron = mcpCollection[0]

        #get energy, momentum and fill truth histos
        trueE = trueNeutron.getEnergy()
        h_truth_E.Fill(trueE)
        dp3 = trueNeutron.getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)
        h_truth_theta.Fill(tlv.Theta())

        #truth entries in photon tree
        E_truth[0] = trueE
        phi_truth[0] = tlv.Phi()
        theta_truth[0] = tlv.Theta()

        #match neutron PFO
        maxpT = -1.
        PFOCollection = event.getCollection('PandoraPFOs')
        for PFO in PFOCollection:
            E_pf = PFO.getEnergy()
            #require minimum PF energy
            if E_pf < 10.:
                continue
            p3 = PFO.getMomentum()
            ptlv = TLorentzVector()
            ptlv.SetPxPyPzE(p3[0], p3[1], p3[2], E_pf)
            pT = ptlv.Perp()
            if pT < maxpT:
                continue
            if PFO.getType() != 2112:
                continue
            theta_temp = ptlv.Theta()
            phi_temp = ptlv.Phi()
            DR = find_dr(phi_temp, theta_temp, phi_truth[0], theta_truth[0])
            if DR > 0.4:
                continue
            else:
                maxpT = pT
                E[0] = E_pf
                theta[0] = ptlv.Theta()
                phi[0] = ptlv.Phi()
             
        if maxpT < 0:
            #no match
            h_reco_neutron_E.Fill(-1.)
            E[0] = -1.
            theta[0] = -1.
            phi[0] = -10.
        else:
            n_pass+=1
            h_reco_neutron_E.Fill(E[0])
            if theta_truth[0] > 0.175 and theta_truth[0] < 2.96:
                h_E_pass.Fill(E_truth[0])
                E_pass_test.append(E_truth[0])
            if E_truth[0] > 10:
                h_theta_pass.Fill(theta_truth[0])
        if theta_truth[0] > 0.175 and theta_truth[0] < 2.96:
            h_E_all.Fill(E_truth[0])
            E_all_test.append(E_truth[0])
        if E_truth[0] > 10:
            h_theta_all.Fill(theta_truth[0])
        neutron_tree.Fill()
        n_all+=1
    reader.close()
# write histograms
output_file = TFile(options.outFile, 'RECREATE')
for histo in histos_list:
    histo.Write()
neutron_tree.Write()
output_file.Close()


