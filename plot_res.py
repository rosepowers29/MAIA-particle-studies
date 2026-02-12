import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import pandas as pd
from argparse import ArgumentParser

# argparse: import (1) No BIB file, if present (2) BIB file, if present (3) particle (4) region (5) outfile/label
parser= ArgumentParser()
parser.add_argument("-n", "--noBIBfile", dest = 'noBIBfile', default=None, help="Path to no-BIB csv")
parser.add_argument("-b", "--BIBfile", dest = 'BIBfile', default=None, help="Path to BIB csv")
parser.add_argument("-p", "--particle", dest = 'particle', default = "Photon", help="particle")
parser.add_argument("-r", "--region", dest = 'region', default = "e", help = "region of detector")
parser.add_argument("-o", "--outfile", dest = 'outfile', default = "E_reso_MAIA_photon_v8_e.pdf", help = "label (version)")
args = parser.parse_args()

def get_arr_and_err(csv_in, arrkey):
    """
    Docstring for get_res
    
    :param csv_in: Description
    """
    df_ = pd.read_csv(csv_in)
    res_arr = df_[arrkey]
    res_err_arr = df_[arrkey+"_err"]
    return res_arr, res_err_arr


def make_labels(ax, region, particle):
    """
    Docstring for make_labels
        
    :param ax: Description
    :param region: Description
    :param particle: Description
    """
    # make region label
    if region == 't':
        ax.get_figure().text(0.165,0.7,"Transition Region ($0.7<\\theta<0.99$ or $2.15<\\theta<2.44$)",fontsize=9)
        ax.set_ylim(0., 1.)
    elif region == 'b':
        ax.get_figure().text(0.165, 0.7, "Central Barrel Region ($0.99<\\theta<2.15$)",fontsize=9)
        ax.set_ylim(0., 1.)
    elif region == 'e':
        ax.get_figure().text(0.165,0.7, "Endcap Region ($0.175<\\theta < 0.7$ or $2.44< \\theta < 2.96$)",fontsize=9)
        ax.set_ylim(0., 1.)

    # make standard MAIA label
    hep.cms.label(exp = "Muon Collider", data = False, label = 'with BIB+IPP (EU24 Lattice) ',
                  rlabel='$MAIA$ Detector Concept', loc=2, italic=(1,0,0), pad=(0.0), ax = ax)
    ax.get_figure().text(0.165,0.74,"$\sqrt{s}=10$ TeV",fontsize=9)

    # make axis labels
    ax.set_xlabel(f"True {particle} Energy [GeV]", loc = 'right')
    ax.set_ylabel(f"{particle} Energy Resolution ($\\sigma /( \\mu$+1)) [%]", loc='top')

def format_plots(ax, big_ticks_x, small_ticks_x, big_ticks_y, small_ticks_y, log=False):
    """
    Docstring for format_plots
     
    :param ax: Description
    :param big_ticks_x: Description
    :param small_ticks_x: Description
    :param big_ticks_y: Description
    :param small_ticks_y: Description
    :param log: Description
    """
    ax.xaxis.set_major_locator(MultipleLocator(big_ticks_x))
    ax.xaxis.set_major_formatter('{x:.1f}')
    ax.xaxis.set_minor_locator(MultipleLocator(small_ticks_x))
     
    ax.yaxis.set_major_locator(MultipleLocator(big_ticks_y))
    ax.yaxis.set_major_formatter('{x:.2f}')
    ax.yaxis.set_minor_locator(MultipleLocator(small_ticks_y))

    ax.grid(axis = 'y', linestyle = ':')
    ax.legend(loc='upper right',frameon=False, fontsize=9)
    if log:
       ax.set_xscale('log')

def plot_data(ax, E_arr, E_err, res_arr, res_err, B=False):
    """
    Docstring for plot_data
     
    :param ax: Description
    :param E_arr: Description
    :param E_err: Description
    :param res_arr: Description
    :param res_err: Description
    :param B: Description
    """
    if B:
        color = "red"
        label = "With BIB+IPP"
    else:
        color = "blue"
        label = "Without BIB+IPP"

    ax.errorbar(E_arr, res_arr, xerr = E_err, yerr = res_err, fmt = ".", color = color, label = label)

def main(infile_nb, infile_b, region, particle, outfile):
    fig, ax = plt.subplots()
    if infile_nb is not None:
        # get files
        E_arr, E_err_arr = get_arr_and_err(infile_nb, "E")
        res_arr, res_err_arr = get_arr_and_err(infile_nb, "res")

        # plot errorbars
        plot_data(ax, E_arr, E_err_arr, res_arr, res_err_arr, False)
    
    if infile_b is not None:
        # get files
        E_arr_B, E_err_arr_B = get_arr_and_err(infile_b, "E")
        res_arr_B, res_err_arr_B = get_arr_and_err(infile_b, "res")
        
        # plot errorbars
        plot_data(ax, E_arr_B, E_err_arr_B, res_arr_B, res_err_arr_B, True)
    
    make_labels(ax, region, particle)
    format_plots(ax, 10, 5, 0.2, 0.05, False)

    fig.savefig(outfile)
    print(outfile,"created")
    plt.close(fig)

if __name__ == "__main__":
    main(args.noBIBfile, args.BIBfile, args.region, args.particle, args.outfile)

