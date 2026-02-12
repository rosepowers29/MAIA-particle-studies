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
parser.add_option('-p', '--particle', default = "photon")
(options, args) = parser.parse_args()

if options.particle == "photon":
    pid = 22
elif options.particle == "neutron":
    pid = 2112
else:
    raise Exception


arrBins_theta = np.linspace(0, TMath.Pi(), 50)
arrBins_E = array('d', (0., 5., 10., 15., 20., 25., 50., 100., 250., 500., 1000.))#, 2500., 5000.))
arrBins_E_lowrange = array('d', (0., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 75., 100., 125., 150., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 900., 1000.))
arrBins_phi = np.linspace(-TMath.Pi(), TMath.Pi(), 50)

# declare histograms

h_PFO_energy = TH1D('PFO_energy', 'PFO_energy', len(arrBins_E_lowrange)-1, arrBins_E_lowrange)
h_PFO_theta = TH1D('PFO_theta', 'PFO_theta', len(arrBins_theta)-1, arrBins_theta)
h_PFO_phi = TH1D('PFO_phi', 'PFO_phi', len(arrBins_phi)-1, arrBins_phi)


# Histo list for writing to outputs
histos_list = [
        h_PFO_energy,
        h_PFO_theta,
        h_PFO_phi
               ]

for histo in histos_list:
    histo.SetDirectory(0)



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

        # get PFO collection ,fill histograms
        PFOCollection = event.getCollection('PandoraPFOs')
        for PFO in PFOCollection:

            if PFO.getType() != pid:
                continue

            E_pf = PFO.getEnergy()
            p3 = PFO.getMomentum() 
            ptlv = TLorentzVector() 
            ptlv.SetPxPyPzE(p3[0], p3[1], p3[2], E_pf) 
            pT = ptlv.Perp() 
            theta = ptlv.Theta() 
            phi = ptlv.Phi() 
            h_PFO_energy.Fill(E_pf)
            h_PFO_theta.Fill(theta)
            h_PFO_phi.Fill(phi)
    
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

output_file.Close()
