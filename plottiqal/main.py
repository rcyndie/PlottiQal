import os
from optparse import OptionParser
import numpy as np
from daskms.experimental.zarr import xds_from_zarr
from plotting_gains import plot_freq_vs_time
from plotting_params import plot_params

#************************************************************************************************************

# COMMAND LINE OPTIONS
def create_parser():
    parser = OptionParser(usage='%prog [options]')
    parser.add_option("--plotop", help="Plotting options: gains, params or all", default="all")
    parser.add_option("--sol1", help="QuartiCal solutions (xarray Dataset)")
    parser.add_option("--sol2", help="Difference with solutions sol", default=None)
    parser.add_option("--ref_ant", help="Reference antenna", default=31)
    parser.add_option("--corr", help="Specify correlation axis", default=0)
    parser.add_option("--scann", help="Specify scan number", default=0)
    parser.add_option("--label0", help="Figure label", default="K")
    parser.add_option("--plotpar", help="Specify parameter to be plotted (delay or TEC)", default="delay")
    parser.add_option("--par_ind", help="Specify parameter index", default=1)
    parser.add_option("--ref_nparr", help="Use to compare parameter plots (same time axis)", default=None)
    parser.add_option("--ftick0", help="Specify frequency axis tick 1", default=50)
    parser.add_option("--ftick1", help="Specify frequency axis tick 2", default=60)
    parser.add_option("--zoomed", help="Additional zoomed in plot on specific colourbar", default=False)
    parser.add_option("-o", "--opdir", help="Output folder to store plots", default="")

    return parser


def main():
    (options,args) = create_parser().parse_args()

    arr0 = xds_from_zarr(options.sol1)
    ref_ant = int(options.ref_ant)
    corr = int(options.corr)
    scann = int(options.scann)
    label0 = options.label0
    plotpar = options.plotpar
    par_ind = int(options.par_ind)
    ftick0 = int(options.ftick0)
    ftick1 = int(options.ftick1)
    zoomed = options.zoomed

    path_plotting = options.opdir + "plots/"
    if not os.path.exists(path_plotting):
        os.makedirs(path_plotting)
        print(f"Directory {path_plotting} created.")
    else:
        print(f"Directory {path_plotting} already exists!")
    

    if options.plotop in ["all", "gains"]:
        #Plotting amplitude and phase.
        f_dict0 = {"func": np.angle, "label1": "ph"}
        f_dict1 = {"func": np.abs, "label1": "amp"}

        #Make the amplitude and phase plots for all terms in the chain.
        for dict1 in [f_dict0, f_dict1]:
            plot_freq_vs_time(path_plotting, arr0, dict1, scann=scann, ref_ant=ref_ant, label0=label0, corr=corr, ftick0=ftick0, \
                ftick1=ftick1, zoomed=False, arr2=None)
            if zoomed:
                plot_freq_vs_time(path_plotting, arr0, dict1, scann=scann, ref_ant=ref_ant, label0=label0, corr=corr, ftick0=ftick0, \
                    ftick1=ftick1, zoomed=options.zoomed, arr2=None)

        if options.sol2 is not None:
            arr2 = xds_from_zarr(options.sol2)
            for dict1 in [f_dict0, f_dict1]:
                plot_freq_vs_time(path_plotting, arr0, dict1, scann=scann, ref_ant=ref_ant, label0="diff_"+label0, corr=corr, ftick0=ftick0, \
                    ftick1=ftick1, zoomed=False, arr2=arr2)
                if zoomed:
                    plot_freq_vs_time(path_plotting, arr0, dict1, scann=scann, ref_ant=ref_ant, label0="diff_"+label0, corr=corr, ftick0=ftick0, \
                        ftick1=ftick1, zoomed=options.zoomed, arr2=arr2)

    if options.plotop in ["all", "params"]:
        plot_params(path_plotting, arr0, scann, ref_ant=ref_ant, label0=label0, plotpar=plotpar, par_ind=par_ind, \
            ref_nparr=None)
        if options.ref_nparr is not None:
            plot_params(path_plotting, arr0, scann, ref_ant=ref_ant, label0=label0, plotpar=plotpar, par_ind=par_ind, \
                ref_nparr=np.load(options.ref_nparr))


#**********************************************************************************************************************
if __name__ == "__main__":
    main()

