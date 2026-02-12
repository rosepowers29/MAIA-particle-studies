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
parser.add_option('-s', '--start', help='index of first file to read',
                type = int, default = '0')
parser.add_option('-e', '--end', help = 'index of last file to read',
                type = int, default = '1000')
(options, args) = parser.parse_args()

arrBins_theta = np.linspace(0.175, TMath.Pi()-0.175, 30)#array('d', (0., 30.*TMath.Pi()/180., 40.*TMath.Pi()/180., 50.*TMath.Pi()/180., 60.*TMath.Pi()/180., 70.*TMath.Pi()/180.,
                            #90.*TMath.Pi()/180., 110.*TMath.Pi()/180., 120.*TMath.Pi()/180., 130.*TMath.Pi()/180., 140.*TMath.Pi()/180., 150.*TMath.Pi()/180., TMath.Pi()))
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000.))#, 2500., 5000.))
arrBins_E_lowrange = array('d', (0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 75., 100., 125., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.))
arrBins_fakes = np.linspace(0,1000,250)

#Global Flag
do_calo_tree = False

# declare histograms

h_E_pass = TH1D('E_pass', 'E_pass', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_E_all = TH1D('E_all', 'E_all', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)

h_theta_pass = TH1D('theta_pass', 'theta_pass', len(arrBins_theta)-1, arrBins_theta)
h_theta_all = TH1D('theta_all', 'theta_all', len(arrBins_theta)-1, arrBins_theta)



# Histo list for writing to outputs
histos_list = [
               h_E_pass, h_E_all,
               h_theta_pass, h_theta_all,
               ]

for histo in histos_list:
    histo.SetDirectory(0)

####################################
photon_tree = TTree("photon_tree", "photon_tree")
E = array('d', [0])
phi = array('d', [0])
theta = array('d', [0])
E_truth = array('d', [0])
phi_truth = array('d', [0])
theta_truth = array('d', [0])
pT = array('d', [0])
pT_truth = array('d', [0])
N_photon = array('d', [0])
photon_tree.Branch("E",  E,  'var/D')
photon_tree.Branch("phi", phi, 'var/D')
photon_tree.Branch("theta", theta, 'var/D')
photon_tree.Branch("E_truth",  E_truth,  'var/D')
photon_tree.Branch("phi_truth", phi_truth, 'var/D')
photon_tree.Branch("theta_truth", theta_truth, 'var/D')
photon_tree.Branch("N_photon", N_photon, 'var/D')


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
    if filenum < options.start:
        filenum = filenum + 1
        continue
    if filenum > options.end:
        break
    # create a reader and open an LCIO file
    reader = IOIMPL.LCFactory.getInstance().createLCReader()
    try:
        reader.open(file)
    except Exception:
        #let it skip the bad files without breaking
        print("skipping file", filenum)
        filenum = filenum+1
        continue
    filenum=filenum+1
    # loop over all events in the file
    for ievt, event in enumerate(reader):
        if filenum % 10 == 0:
            print(" ")
            print("File "+str(filenum))
            print("Processing event " + str(ievt))

        #print(event.getCollectionNames())

        # Fill the truth-level histos, the first particle is always the gun
        mcpCollection = event.getCollection('MCParticle')
        truePhoton = mcpCollection[0]
        
        #get energy, momentum and fill truth histos
        trueE = truePhoton.getEnergy()
        dp3 = truePhoton.getMomentum()
        tlv = TLorentzVector()
        tlv.SetPxPyPzE(dp3[0], dp3[1], dp3[2], trueE)

        #truth entries in photon tree
        E_truth[0] = trueE
        phi_truth[0] = tlv.Phi()
        theta_truth[0] = tlv.Theta()

        #match photon PFO
        maxpT = -1.
        PFOCollection = event.getCollection('PandoraPFOs')
        for PFO in PFOCollection:
            E_pf = PFO.getEnergy()
            p3 = PFO.getMomentum()
            ptlv = TLorentzVector()
            ptlv.SetPxPyPzE(p3[0], p3[1], p3[2], E_pf)
            pT = ptlv.Perp()
            if pT < maxpT:
                #print("failed: pT too small")
                continue
            if PFO.getType() != 22:
                #print("failed: wrong type")
                continue
            theta_temp = ptlv.Theta()
            phi_temp = ptlv.Phi()
            DR = find_dr(phi_temp, theta_temp, phi_truth[0], theta_truth[0])
            if DR > 0.1:
                #print("failed: dR>0.1")
                continue
            else:
                maxpT = pT
                E[0] = E_pf
                theta[0] = ptlv.Theta()
                phi[0] = ptlv.Phi()
        
        if maxpT < 0:
            #no match
            E[0] = -1.
            theta[0] = -1.
            phi[0] = -10. 
        else:
            if theta_truth[0] > 0.175 and theta_truth[0] < 2.96:
                h_E_pass.Fill(E_truth[0])
            if E_truth[0] > 10:
                h_theta_pass.Fill(theta_truth[0])
        if theta_truth[0] > 0.175 and theta_truth[0] < 2.96:
            h_E_all.Fill(E_truth[0])
        if E_truth[0] > 10:
            h_theta_all.Fill(theta_truth[0])
        photon_tree.Fill()
    reader.close()
    
# write histograms
print("out of main loop")
try:
    outfile_name = f"{options.outFile}_{options.start}-{options.end}.root"
except Exception:
    outfile_name = "outfile_FIXME.root"
output_file = TFile(outfile_name, 'RECREATE')
for histo in histos_list:
    histo.Write()
photon_tree.Write()

output_file.Close()
